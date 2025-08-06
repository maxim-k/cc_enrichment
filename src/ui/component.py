import json

import streamlit as st


def d3_force_graph(json_obj, height=500):
    """
    Render a D3.js force-directed graph inside Streamlit given a JSON object.
    Expects format:
    {
      "nodes": [
         {"id": "...", "type": "term" or "gene", "label": "...", "color": "#hex"}  # term has color
      ],
      "links": [
         {"source": "...", "target": "...", "color": "#hex"}
      ]
    }
    """
    data_json = json.dumps(json_obj).replace("</", "<\/")
    html = f"""
    <!DOCTYPE html>
    <html>
    <head>
      <meta charset="utf-8">
      <style>
        body {{ margin:0; font-family: system-ui, -apple-system, BlinkMacSystemFont, sans-serif; }}
        #chart {{ width:100%; height:100%; position: relative; }}
        .tooltip {{
          position: absolute;
          padding: 6px 10px;
          background: rgba(0,0,0,0.75);
          color: white;
          font-size: 12px;
          border-radius: 4px;
          pointer-events: none;
          z-index: 10;
          white-space: nowrap;
        }}
        .legend-item {{
          display: flex;
          align-items: center;
          gap: 6px;
          font-size: 12px;
          margin-right: 12px;
        }}
        .legend-swatch {{
          width: 12px;
          height: 12px;
          border-radius: 2px;
          display: inline-block;
        }}
      </style>
      <script src="https://d3js.org/d3.v7.min.js"></script>
    </head>
    <body>
      <div id="chart"></div>
      <script>
        const graph = {data_json};

        const container = d3.select("#chart");
        const width = container.node().clientWidth || 800;
        const height = {height};

        // tooltip
        const tooltip = d3.select("body")
          .append("div")
          .attr("class", "tooltip")
          .style("opacity", 0);

        const svg = container
          .append("svg")
          .attr("width", width)
          .attr("height", height);

        // define simulation
        const simulation = d3.forceSimulation(graph.nodes)
          .force("link", d3.forceLink(graph.links).id(d => d.id).distance(80).strength(0.8))
          .force("charge", d3.forceManyBody().strength(-120))
          .force("center", d3.forceCenter(width / 2, height / 2))
          .force("collision", d3.forceCollide().radius(d => (d.type === "term" ? 14 : 8) + 4));

        // run a few iterations to stabilize initial layout
        for (let i = 0; i < 120; ++i) simulation.tick();

        // draw links
        const link = svg.append("g")
          .attr("stroke-width", 1.5)
          .selectAll("line")
          .data(graph.links)
          .join("line")
            .attr("stroke", d => d.color || "#aaa")
            .attr("stroke-opacity", 0.8);

        // draw nodes
        const node = svg.append("g")
          .selectAll("circle")
          .data(graph.nodes)
          .join("circle")
            .attr("r", d => d.type === "term" ? 12 : 6)
            .attr("fill", d => d.type === "term" ? (d.color || "#444") : "#999")
            .attr("stroke", "#222")
            .attr("stroke-width", d => d.type === "term" ? 1.5 : 1)
            .call(drag(simulation));

        // labels for term nodes only
        const labels = svg.append("g")
          .selectAll("text")
          .data(graph.nodes)
          .join("text")
            .text(d => d.label || d.id)
            .attr("font-family", "Arial")
            .attr("font-size", 14)
            .attr("font-weight", "bold")
            .attr("stroke", "#FFFFFF80")

            .attr("dy", d => d.type === "term" ? "-1em" : "-0.5em")
            .attr("text-anchor", "middle")
            .attr("pointer-events", "none");

        // hover interactions
        node.on("mouseover", (event, d) => {{
            tooltip
              .style("opacity", 1)
              .html(`<strong>${{d.label || d.id}}</strong><br/>type: ${{d.type}}`)
              .style("left", (event.pageX + 8) + "px")
              .style("top", (event.pageY + 8) + "px");
        }}).on("mousemove", (event) => {{
            tooltip
              .style("left", (event.pageX + 8) + "px")
              .style("top", (event.pageY + 8) + "px");
        }}).on("mouseout", () => {{
            tooltip.style("opacity", 0);
        }});

        // simulation tick update
        simulation.on("tick", () => {{
          link
            .attr("x1", d => d.source.x)
            .attr("y1", d => d.source.y)
            .attr("x2", d => d.target.x)
            .attr("y2", d => d.target.y);

          node
            .attr("cx", d => d.x)
            .attr("cy", d => d.y);

          labels
            .attr("x", d => d.x)
            .attr("y", d => d.y);
        }});

        // drag helpers
        function drag(sim) {{
          function dragstarted(event, d) {{
            if (!event.active) sim.alphaTarget(0.3).restart();
            d.fx = d.x;
            d.fy = d.y;
          }}

          function dragged(event, d) {{
            d.fx = event.x;
            d.fy = event.y;
          }}

          function dragended(event, d) {{
            if (!event.active) sim.alphaTarget(0);
            d.fx = null;
            d.fy = null;
          }}

          return d3.drag()
            .on("start", dragstarted)
            .on("drag", dragged)
            .on("end", dragended);
        }}

        // legend for terms
        const termColors = Array.from(new Set(graph.nodes
          .filter(n => n.type === "term" && n.color)
          .map(n => n.color)));
        if (termColors.length) {{
          const legend = container.append("div")
            .style("position", "absolute")
            .style("top", "8px")
            .style("left", "8px")
            .style("display", "flex");
          termColors.forEach(c => {{
            const item = legend.append("div").attr("class", "legend-item");
            item.append("div")
              .attr("class", "legend-swatch")
              .style("background", c);
            item.append("div").text("term " + c);
          }});
        }}

      </script>
    </body>
    </html>
    """
    st.components.v1.html(html, height=height, scrolling=True)
