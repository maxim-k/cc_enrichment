from typing import Dict, List

import networkx as nx
import plotly.graph_objects as go
import pydot
from plotly.colors import sample_colorscale


def _sample_scale_hex(scale_name: str, n: int) -> List[str]:
    """Evenly sample a Plotly continuous colorscale and return hex colors.

    Parameters
    ----------
    scale_name : str
        Plotly continuous colorscale name, e.g. "Viridis", "Turbo", "Plasma".
    n : int
        Number of colors to sample.

    Returns
    -------
    List[str]
        Colors like '#1f77b4'.
    """
    if n <= 0:
        return []
    if n == 1:
        stops = [0.5]
    else:
        stops = [i / (n - 1) for i in range(n)]

    rgb_list = sample_colorscale(scale_name, stops)

    def rgb_to_hex(rgb: str) -> str:
        # rgb is like 'rgb(12, 34, 56)'
        rgb_vals = rgb.strip()[4:-1].split(",")  # ['12', ' 34', ' 56']
        r, g, b = (int(v) for v in rgb_vals)
        return f"#{r:02x}{g:02x}{b:02x}"

    return [rgb_to_hex(rgb) for rgb in rgb_list]


def _colorize_term_nodes(dot_src: str, hex_color: str) -> str:
    """Insert fillcolor for *term* nodes using simple string ops, no regex parsing.

    Assumes term lines follow the fixed format:
    "term_..." [label="...", style=filled, fontcolor="white"];

    We replace the substring 'style=filled,' with
    'style=filled, fillcolor="<hex_color>",' only for lines that start with a term node.
    """
    out_lines: List[str] = []

    for line in dot_src.splitlines():
        if 'type="term"' in line and 'fillcolor=' not in line:
            if 'style=filled,' in line:
                line = line.replace(
                    'style=filled,', f'style=filled, fillcolor="{hex_color}",'
                )
            elif 'style=filled' in line:
                line = line.replace(
                    'style=filled', f'style=filled, fillcolor="{hex_color}"'
                )
            elif 'fontcolor=' in line:
                line = line.replace(
                    'fontcolor=', f'fillcolor="{hex_color}", fontcolor='
                )
            else:
                idx = line.rfind("]")
                if idx != -1:
                    line = line[:idx] + f', fillcolor="{hex_color}"' + line[idx:]
        out_lines.append(line)
    return "\n".join(out_lines)


def _build_library_legend(lib_to_color: Dict[str, str]) -> str:
    """Create a legend subgraph mapping each library to its color."""
    lines = [
        "  subgraph cluster_legend_libs {",
        '    label="Libraries";',
        '    node [shape=box, style=filled, fontcolor="white"];',
    ]
    for i, (lib, color) in enumerate(lib_to_color.items(), 1):
        safe_id = f"legend_lib_{i}"
        lbl = lib.replace('"', '\\"')
        lines.append(f'    {safe_id} [label="{lbl}", fillcolor="{color}"];')
    lines.append("  }")
    return "\n".join(lines)


def merge_iterative_dot(
    per_lib_dots: Dict[str, str], scale_name: str = "Viridis"
) -> str:
    """Merge per-library DOT strings and color term nodes per library.

    per_lib_dots : {library_name: dot_string}
    scale_name   : Plotly continuous scale name.
    """
    libs = list(per_lib_dots.keys())
    colors = _sample_scale_hex(scale_name, len(libs))
    lib_to_color = dict(zip(libs, colors))

    # 1) Colorize each library snippet
    colorized = {
        lib: _colorize_term_nodes(dot_src, lib_to_color[lib])
        for lib, dot_src in per_lib_dots.items()
    }

    # 2) Strip outer graph braces and collect unique lines
    def strip_outer(snippet: str) -> str:
        # naive: find first '{' and last '}'
        start = snippet.find("{")
        end = snippet.rfind("}")
        if start != -1 and end != -1 and end > start:
            return snippet[start + 1 : end]
        return snippet

    seen = set()
    body_lines: List[str] = []
    graph_attrs: List[str] = []
    node_defaults: str | None = None

    for snip in colorized.values():
        inner = strip_outer(snip)
        for ln in inner.splitlines():
            stripped = ln.strip()
            if not stripped:
                continue
            if stripped.startswith("graph ") and stripped.endswith("];"):
                if stripped not in graph_attrs:
                    graph_attrs.append(stripped)
                continue
            if stripped.startswith("node ") and stripped.endswith("];"):
                if node_defaults is None:
                    node_defaults = stripped
                continue
            if stripped not in seen:
                seen.add(stripped)
                body_lines.append(ln)

    # 3) Remove any old cluster_legend_libs block
    cleaned: List[str] = []
    skip = False
    for ln in body_lines:
        if "subgraph cluster_legend_libs" in ln:
            skip = True
        if skip and "}" in ln:
            skip = False
            continue
        if not skip:
            cleaned.append(ln)
    body_lines = cleaned

    # 4) Assemble final DOT
    out = ["graph iterative_enrichment_all {"]
    for l in graph_attrs:
        out.append(f"  {l}")
    if node_defaults:
        out.append(f"  {node_defaults}")

    out.append(_build_library_legend(lib_to_color))

    for l in body_lines:
        out.append(f"  {l}")

    out.append("}")
    return "\n".join(out)


def dot_to_plotly(
    dot_input: str,
    node_size: int = 7,
    edge_width: int = 1,
    layout_k: float = 0.7,
    layout_iterations: int = 100,
) -> go.Figure:
    """
    Convert a DOT file into an interactive Plotly network figure.

    Parameters:
    - dot_input: path to the .dot file
    - node_size: base marker size for nodes
    - edge_width: line width for edges
    - layout_k: optimal distance between nodes for spring layout
    - layout_iterations: iterations for force-directed layout

    Returns:
    - fig: plotly.graph_objects.Figure
    """
    # 1. Parse DOT with pydot
    try:
        graphs = pydot.graph_from_dot_data(dot_input)
    except pydot.PydotException as e:
        raise ValueError(f"Failed to parse DOT data: {dot_input!r}") from e
    if not graphs:
        raise ValueError(f"Failed to parse DOT data: {dot_input!r}")

    dot = graphs[0]

    # 2. Convert to NetworkX graph
    G = nx.Graph()
    for node in dot.get_nodes():
        name = node.get_name().strip('"')
        # skip meta or empty nodes
        if not name or name.lower() == "node" or name.lower() == "graph":
            continue
        attrs = node.get_attributes() or {}
        label = attrs.get("label", name).strip('"')
        raw = attrs.get("fillcolor")
        color = raw.strip('"') if raw else "#888888"
        G.add_node(name, label=label, color=color)
    for edge in dot.get_edges():
        src = edge.get_source().strip('"')
        dst = edge.get_destination().strip('"')
        if src and dst and G.has_node(src) and G.has_node(dst):
            G.add_edge(src, dst)

    # 3. Compute positions with force-directed layout
    pos = nx.spring_layout(
        G,
        k=layout_k,
        iterations=layout_iterations,
    )

    # 4. Build edge trace
    edge_x, edge_y = [], []
    for u, v in G.edges():
        x0, y0 = pos[u]
        x1, y1 = pos[v]
        edge_x += [x0, x1, None]
        edge_y += [y0, y1, None]
    edge_trace = go.Scatter(
        x=edge_x,
        y=edge_y,
        mode="lines",
        line=dict(width=edge_width, color="rgba(150,150,150,0.3)"),
        hoverinfo="none",
    )

    # 5. Build node trace
    node_x, node_y, node_text, node_color = [], [], [], []
    for n, data in G.nodes(data=True):
        x, y = pos[n]
        node_x.append(x)
        node_y.append(y)
        node_text.append(data.get("label", n))
        # use fillcolor if provided, else default
        fill = data.get("color")
        node_color.append(fill if fill else "rgba(50,50,250,0.6)")

    node_trace = go.Scatter(
        x=node_x,
        y=node_y,
        mode="markers+text",
        text=node_text,
        textposition="top center",
        marker=dict(
            size=node_size,
            color=node_color,
            opacity=0.7,
            line=dict(width=1, color="rgba(0,0,0,0.2)"),
        ),
        hoverinfo="text",
    )

    # 6. Assemble figure
    fig = go.Figure(data=[edge_trace, node_trace])
    fig.update_layout(
        margin=dict(l=0, r=0, t=0, b=0),
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        plot_bgcolor="white",
        autosize=False,
        width=1000,
        height=1000,
    )
    return fig
