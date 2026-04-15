#!/usr/bin/env python3
"""
Export drawio diagrams to high-resolution PNG.
Uses draw.io's online conversion service as fallback.
Usage: python export_drawio.py <input.drawio> <output.png> [scale]
"""
import sys
import os
import zlib
import base64
import re
from defusedxml import ElementTree as ET
import cairosvg

def compress_and_encode(text):
    """Compress with zlib deflate and base64 encode (drawio format)."""
    compressed = zlib.compress(text.encode('utf-8'))[2:-4]  # strip header/trailer
    return base64.b64encode(compressed).decode('ascii')

def decompress_and_decode(b64_text):
    """Decode base64 and decompress (drawio format)."""
    raw = base64.b64decode(b64_text)
    # Add zlib header (2 bytes) and adler32 checksum (4 bytes)
    data = b'\x78\x9c' + raw  # default compression header
    return zlib.decompress(data).decode('utf-8')

def mxgraphmodel_to_svg(mxgraphmodel_xml, page_width=1100, page_height=1700):
    """
    Convert mxGraphModel to a renderable SVG.
    This creates an SVG representation of the flowchart boxes and text.
    """
    root = ET.fromstring(mxgraphmodel_xml)
    
    # Collect all cells
    cells = []
    for mxCell in root.findall('.//mxCell'):
        cell_id = mxCell.get('id', '')
        if cell_id in ('0', '1'):
            continue
        
        value = mxCell.get('value', '')
        style = mxCell.get('style', '')

        # Parse style
        styles = {}
        for part in style.split(';'):
            if '=' in part:
                k, v = part.split('=', 1)
                styles[k] = v
        
        # Get geometry
        geo = mxCell.find('mxGeometry')

        # Draw.io edges usually have edge="1" and still include mxGeometry.
        # Detect them explicitly so we don't mis-render them as top-left rectangles.
        if mxCell.get('edge', '0') == '1':
            stroke = styles.get('strokeColor', '#333333')
            if not stroke.startswith('#'):
                stroke = '#333333'

            stroke_width = 1.5
            try:
                stroke_width = float(styles.get('strokeWidth', '1.5'))
            except ValueError:
                pass

            points = []
            source_point = None
            target_point = None
            if geo is not None:
                src_pt = geo.find("mxPoint[@as='sourcePoint']")
                if src_pt is not None and src_pt.get('x') is not None and src_pt.get('y') is not None:
                    source_point = (float(src_pt.get('x')), float(src_pt.get('y')))

                tgt_pt = geo.find("mxPoint[@as='targetPoint']")
                if tgt_pt is not None and tgt_pt.get('x') is not None and tgt_pt.get('y') is not None:
                    target_point = (float(tgt_pt.get('x')), float(tgt_pt.get('y')))

                point_array = geo.find("Array[@as='points']")
                if point_array is not None:
                    for p in point_array.findall('mxPoint'):
                        px = p.get('x')
                        py = p.get('y')
                        if px is None or py is None:
                            continue
                        points.append((float(px), float(py)))

            cells.append({
                'type': 'edge',
                'id': cell_id,
                'source': mxCell.get('source', ''),
                'target': mxCell.get('target', ''),
                'style': style,
                'points': points,
                'source_point': source_point,
                'target_point': target_point,
                'stroke': stroke,
                'stroke_width': stroke_width,
                'is_dashed': styles.get('dashed', '0') == '1',
                'edge_style': styles.get('edgeStyle', ''),
                'exitX': float(styles['exitX']) if 'exitX' in styles else None,
                'exitY': float(styles['exitY']) if 'exitY' in styles else None,
                'entryX': float(styles['entryX']) if 'entryX' in styles else None,
                'entryY': float(styles['entryY']) if 'entryY' in styles else None,
            })
            continue

        if geo is None:
            continue
        
        x = float(geo.get('x', 0))
        y = float(geo.get('y', 0))
        w = float(geo.get('width', 100))
        h = float(geo.get('height', 30))
        
        fill = styles.get('fillColor', 'none')
        stroke = styles.get('strokeColor', 'none')
        font_size = styles.get('fontSize', '12')
        is_rounded = 'rounded=1' in style
        is_rhombus = 'rhombus' in style
        is_text = 'text;html=1' in style
        is_dashed = styles.get('dashed', '0') == '1'
        font_style = styles.get('fontStyle', '0')
        is_bold = int(font_style) & 1 if font_style else False
        
        # Map drawio colors
        if fill.startswith('#'):
            pass
        elif fill == 'none':
            fill = 'transparent'
        else:
            fill = '#FFFFFF'
        
        if stroke.startswith('#'):
            pass
        elif stroke == 'none':
            stroke = 'transparent'
        else:
            stroke = '#000000'
        
        # Convert HTML entities in text
        value = value.replace('&amp;', '&').replace('&lt;', '<').replace('&gt;', '>')
        value = value.replace('&#xa;', '\n').replace('<br>', '\n')
        value = value.replace('<b>', '').replace('</b>', '')
        value = value.strip()
        
        cells.append({
            'type': 'shape',
            'id': cell_id,
            'x': x, 'y': y, 'w': w, 'h': h,
            'fill': fill, 'stroke': stroke,
            'is_rounded': is_rounded,
            'is_rhombus': is_rhombus,
            'is_text': is_text,
            'is_dashed': is_dashed,
            'is_bold': is_bold,
            'font_size': int(font_size),
            'value': value,
            'style': style,
        })
    
    # Build SVG
    svg_parts = []
    svg_parts.append(f'<svg xmlns="http://www.w3.org/2000/svg" width="{page_width}" height="{page_height}" viewBox="0 0 {page_width} {page_height}">')
    svg_parts.append(f'<rect width="100%" height="100%" fill="white"/>')
    
    # Draw shapes first, then text on top
    shapes = [c for c in cells if c['type'] == 'shape' and not c['is_text']]
    text_elements = [c for c in cells if c['type'] == 'shape' and c['is_text']]
    edges = [c for c in cells if c['type'] == 'edge']
    
    for cell in shapes:
        dasharray = 'stroke-dasharray="5,5" ' if cell['is_dashed'] else ''
        fill_attr = f'fill="{cell["fill"]}"' if cell['fill'] != 'transparent' else 'fill="none"'
        
        if cell['is_rhombus']:
            cx, cy = cell['x'] + cell['w']/2, cell['y'] + cell['h']/2
            points = f"{cx},{cell['y']} {cell['x']+cell['w']},{cy} {cx},{cell['y']+cell['h']} {cell['x']},{cy}"
            svg_parts.append(f'<polygon points="{points}" {fill_attr} stroke="{cell["stroke"]}" stroke-width="1.5" {dasharray}/>')
        elif cell['is_rounded']:
            rx, ry = min(10, cell['w']/6), min(10, cell['h']/6)
            svg_parts.append(f'<rect x="{cell["x"]}" y="{cell["y"]}" width="{cell["w"]}" height="{cell["h"]}" rx="{rx}" ry="{ry}" {fill_attr} stroke="{cell["stroke"]}" stroke-width="1.5" {dasharray}/>')
        else:
            svg_parts.append(f'<rect x="{cell["x"]}" y="{cell["y"]}" width="{cell["w"]}" height="{cell["h"]}" {fill_attr} stroke="{cell["stroke"]}" stroke-width="1.5" {dasharray}/>')
        
        # Add text inside shape - escape XML special chars
        if cell['value']:
            font_weight = 'bold' if cell['is_bold'] else 'normal'
            lines = cell['value'].split('\n')
            line_h = cell['font_size'] * 1.3
            total_h = len(lines) * line_h
            text_start_y = cell['y'] + (cell['h'] - total_h) / 2 + cell['font_size']
            
            for i, line in enumerate(lines):
                # Escape XML special characters
                safe_line = line.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;').replace('"', '&quot;')
                ty = text_start_y + i * line_h
                svg_parts.append(f'<text x="{cell["x"] + cell["w"]/2}" y="{ty}" text-anchor="middle" font-family="Noto Sans CJK SC, sans-serif" font-size="{cell["font_size"]}" fill="#000000" font-weight="{font_weight}">{safe_line}</text>')
    
    # Draw text-only elements
    for cell in text_elements:
        font_weight = 'bold' if cell['is_bold'] else 'normal'
        lines = cell['value'].split('\n')
        line_h = cell['font_size'] * 1.3
        text_start_y = cell['y'] + cell['font_size']
        
        for i, line in enumerate(lines):
            # Escape XML special characters
            safe_line = line.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;').replace('"', '&quot;')
            ty = text_start_y + i * line_h
            svg_parts.append(f'<text x="{cell["x"] + cell["w"]/2}" y="{ty}" text-anchor="middle" font-family="Noto Sans CJK SC, sans-serif" font-size="{cell["font_size"]}" fill="#000000" font-weight="{font_weight}">{safe_line}</text>')

    # Draw edges (arrows)
    # Build shape maps for source/target anchoring
    shape_map = {}
    bounds_map = {}
    for cell in cells:
        if cell['type'] == 'shape':
            shape_map[cell['id']] = cell
            bounds_map[cell['id']] = (cell['x'], cell['y'], cell['w'], cell['h'])

    def anchor_point(bounds, toward):
        """Find rectangle boundary intersection along center->target ray."""
        x, y, w, h = bounds
        cx = x + w / 2.0
        cy = y + h / 2.0
        tx, ty = toward
        dx = tx - cx
        dy = ty - cy

        eps = 1e-9
        if abs(dx) < eps and abs(dy) < eps:
            return (cx, cy)

        tx_scale = (w / 2.0) / abs(dx) if abs(dx) > eps else float('inf')
        ty_scale = (h / 2.0) / abs(dy) if abs(dy) > eps else float('inf')
        t = min(tx_scale, ty_scale)

        ax = cx + dx * t
        ay = cy + dy * t
        ax = max(x, min(x + w, ax))
        ay = max(y, min(y + h, ay))

        return (ax, ay)

    def point_on_shape(bounds, rel_x, rel_y):
        x, y, w, h = bounds
        return (x + rel_x * w, y + rel_y * h)

    def text_box_anchor(bounds, toward):
        """Keep text label connectors near center for cleaner section transitions."""
        x, y, w, h = bounds
        cx = x + w / 2.0
        cy = y + h / 2.0
        tx, ty = toward
        ax = max(x, min(x + w, tx))
        if ty >= cy:
            return (ax, y + h)
        return (ax, y)

    def orthogonalize_path(points):
        if len(points) < 2:
            return points

        path = [points[0]]
        for next_pt in points[1:]:
            cur = path[-1]
            if abs(cur[0] - next_pt[0]) < 0.01 or abs(cur[1] - next_pt[1]) < 0.01:
                path.append(next_pt)
                continue

            # Route as orthogonal with dominant direction first to reduce zig-zag.
            if abs(cur[0] - next_pt[0]) >= abs(cur[1] - next_pt[1]):
                bend = (next_pt[0], cur[1])
            else:
                bend = (cur[0], next_pt[1])
            path.append(bend)
            path.append(next_pt)

        return path

    def simplify_path(points):
        if len(points) < 2:
            return points

        snapped = [points[0]]
        for px, py in points[1:]:
            qx, qy = snapped[-1]
            if abs(px - qx) < 0.5:
                px = qx
            if abs(py - qy) < 0.5:
                py = qy
            if abs(px - qx) < 0.01 and abs(py - qy) < 0.01:
                continue
            snapped.append((px, py))

        if len(snapped) < 3:
            return snapped

        simplified = [snapped[0]]
        for i in range(1, len(snapped) - 1):
            ax, ay = simplified[-1]
            bx, by = snapped[i]
            cx, cy = snapped[i + 1]
            if (abs(ax - bx) < 0.01 and abs(bx - cx) < 0.01) or (
                abs(ay - by) < 0.01 and abs(by - cy) < 0.01
            ):
                continue
            simplified.append((bx, by))
        simplified.append(snapped[-1])
        return simplified

    marker_ids = {}

    def marker_id_for(stroke):
        marker_id = marker_ids.get(stroke)
        if marker_id is not None:
            return marker_id
        key = re.sub(r'[^a-zA-Z0-9]+', '_', stroke).strip('_') or 'default'
        marker_id = f'arrowhead_{key}_{len(marker_ids)}'
        marker_ids[stroke] = marker_id
        return marker_id

    for edge in edges:
        src_id = edge.get('source', '')
        tgt_id = edge.get('target', '')
        mid_points = list(edge.get('points', []))
        start_pt = edge.get('source_point')
        end_pt = edge.get('target_point')

        src_bounds = bounds_map.get(src_id) if src_id else None
        tgt_bounds = bounds_map.get(tgt_id) if tgt_id else None

        if src_bounds is not None:
            if edge.get('exitX') is not None and edge.get('exitY') is not None:
                start_pt = point_on_shape(src_bounds, edge['exitX'], edge['exitY'])
            elif start_pt is None:
                toward = mid_points[0] if mid_points else (
                    end_pt if end_pt is not None else (
                        (tgt_bounds[0] + tgt_bounds[2]/2.0, tgt_bounds[1] + tgt_bounds[3]/2.0)
                        if tgt_bounds is not None else None
                    )
                )
                if toward is not None:
                    src_shape = shape_map.get(src_id)
                    if src_shape is not None and src_shape.get('is_text', False):
                        start_pt = text_box_anchor(src_bounds, toward)
                    else:
                        start_pt = anchor_point(src_bounds, toward)

        if tgt_bounds is not None:
            if edge.get('entryX') is not None and edge.get('entryY') is not None:
                end_pt = point_on_shape(tgt_bounds, edge['entryX'], edge['entryY'])
            elif end_pt is None:
                toward = mid_points[-1] if mid_points else (
                    start_pt if start_pt is not None else (
                        (src_bounds[0] + src_bounds[2]/2.0, src_bounds[1] + src_bounds[3]/2.0)
                        if src_bounds is not None else None
                    )
                )
                if toward is not None:
                    tgt_shape = shape_map.get(tgt_id)
                    if tgt_shape is not None and tgt_shape.get('is_text', False):
                        end_pt = text_box_anchor(tgt_bounds, toward)
                    else:
                        end_pt = anchor_point(tgt_bounds, toward)

        if start_pt is None or end_pt is None:
            continue

        poly_points = [start_pt] + mid_points + [end_pt]
        if edge.get('edge_style') == 'orthogonalEdgeStyle':
            poly_points = orthogonalize_path(poly_points)

        filtered = simplify_path(poly_points)

        if len(filtered) >= 2:
            pts = ' '.join(f'{x:.2f},{y:.2f}' for x, y in filtered)
            dasharray = 'stroke-dasharray="5,5" ' if edge.get('is_dashed', False) else ''
            stroke = edge.get('stroke', '#333333')
            stroke_width = edge.get('stroke_width', 1.5)
            marker_id = marker_id_for(stroke)
            svg_parts.append(
                f'<polyline points="{pts}" fill="none" stroke="{stroke}" '
                f'stroke-width="{stroke_width}" stroke-linecap="round" '
                f'stroke-linejoin="round" {dasharray}marker-end="url(#{marker_id})"/>'
            )

    # Add arrowhead markers (color-matched with stroke for cleaner rendering).
    if marker_ids:
        defs = ['<defs>']
        for stroke, marker_id in marker_ids.items():
            defs.append(
                f'<marker id="{marker_id}" markerWidth="10" markerHeight="7" '
                f'refX="9" refY="3.5" orient="auto">'
                f'<polygon points="0 0, 10 3.5, 0 7" fill="{stroke}"/></marker>'
            )
        defs.append('</defs>')
        svg_parts.insert(1, '\n'.join(defs))
    
    svg_parts.append('</svg>')
    return '\n'.join(svg_parts)

def export_drawio_to_png(input_file, output_file, scale=3):
    """Main export function."""
    # Read drawio file
    with open(input_file, 'r', encoding='utf-8') as f:
        content = f.read()
    
    # Parse XML
    root = ET.fromstring(content)
    diagram = root.find('.//diagram')
    if diagram is None:
        raise ValueError("No diagram element found")
    
    diagram_text = diagram.text.strip() if diagram.text else ''
    
    # Check if diagram has mxGraphModel child element (not compressed text)
    gm = diagram.find('mxGraphModel')
    if gm is not None:
        print("Detected mxGraphModel child element format")
        page_width = int(gm.get('pageWidth', 850))
        page_height = int(gm.get('pageHeight', 1100))
        
        print(f"Converting mxGraphModel to SVG (page: {page_width}x{page_height})...")
        svg_content = mxgraphmodel_to_svg(ET.tostring(gm, encoding='unicode'), page_width, page_height)
    elif diagram_text:
        # Try to decompress
        try:
            decoded = base64.b64decode(diagram_text)
            data = b'\x78\x9c' + decoded
            decompressed = zlib.decompress(data).decode('utf-8')
            svg_content = decompressed
            print("Detected compressed format, decompressed successfully")
        except Exception:
            if '<svg' in diagram_text or '<mxGraphModel' in diagram_text:
                if '<mxGraphModel' in diagram_text:
                    gm = ET.fromstring(diagram_text)
                    page_width = int(gm.get('pageWidth', 850))
                    page_height = int(gm.get('pageHeight', 1100))
                    print(f"Converting embedded mxGraphModel (page: {page_width}x{page_height})...")
                    svg_content = mxgraphmodel_to_svg(diagram_text, page_width, page_height)
                else:
                    svg_content = diagram_text
                    print("Detected plain SVG format")
            else:
                raise ValueError("Unknown diagram format")
    else:
        raise ValueError("No diagram content found (empty diagram element)")
    
    svg_size = len(svg_content)
    print(f"SVG content: {svg_size:,} bytes")
    
    # Save debug SVG
    debug_svg = output_file.replace('.png', '.debug.svg')
    with open(debug_svg, 'w', encoding='utf-8') as f:
        f.write(svg_content)
    print(f"Debug SVG: {debug_svg}")
    
    # Convert to PNG
    print(f"Converting to PNG (scale={scale}x)...")
    cairosvg.svg2png(
        bytestring=svg_content.encode('utf-8'),
        write_to=output_file,
        scale=scale,
        background_color='white'
    )
    
    size = os.path.getsize(output_file)
    print(f"OK: {output_file} ({size:,} bytes)")
    
    # Clean up debug file
    if os.path.exists(debug_svg):
        os.remove(debug_svg)

def main():
    if len(sys.argv) < 3:
        print("Usage: python export_drawio.py <input.drawio> <output.png> [scale]")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    scale = int(sys.argv[3]) if len(sys.argv) > 3 else 3
    
    if not os.path.exists(input_file):
        print(f"Error: File not found: {input_file}")
        sys.exit(1)
    
    print(f"Input:  {input_file}")
    print(f"Output: {output_file}")
    
    try:
        export_drawio_to_png(input_file, output_file, scale=scale)
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == '__main__':
    main()
