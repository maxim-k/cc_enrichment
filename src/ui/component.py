import json
import streamlit as st


def d3_force_graph(json_obj, width=1000, height=800):
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
    data_json = json.dumps(json_obj).replace("</", "<\\/>")
    html = f"""
    <!DOCTYPE html>
    <html>
    <head>
      <meta charset="utf-8">
      <style>
        body {{ margin:0; font-family: system-ui, -apple-system, BlinkMacSystemFont, sans-serif; }}
        #chart {{
            width:98%;
            height:100%;
            position: relative;
        }}
        .chart-container {{
          display: flex;
          justify-content: flex-end; /* push contents to the right */
          align-items: center;       /* vertically center if there’s height */
        }}

        /* reset the button, give it the opacity behavior */
        .icon-button {{
          background: none;
          border: none;
          padding: 4px;      /* tweak as needed around the SVG */
          margin: 0;
          opacity: 0.5;      /* default half-transparent */
          cursor: pointer;
        }}

        .icon-button:hover,
        .icon-button:focus {{
          opacity: 1;        /* fully opaque on hover/focus */
        }}

        /* optional: ensure SVG fills its container and inherits color */
        .icon-button svg {{
          display: block;       /* kills inline-svg descenders */
          fill: currentColor;   /* svg path will pick up text color */
          height: 1em;
          width: 1em;
        }}
</style>        
      </style>
      <script src="https://d3js.org/d3.v7.min.js"></script>
      <script src="https://cdn.jsdelivr.net/npm/save-svg-as-png@1.4.17/lib/saveSvgAsPng.js"></script>
    </head>
    <body>
        <div id="chart-container">
            <button id="download-png" class="icon-button" aria-label="Download PNG">
                <svg viewBox="0 0 1000 1000" class="icon" height="1em" width="1em"><path d="m500 450c-83 0-150-67-150-150 0-83 67-150 150-150 83 0 150 67 150 150 0 83-67 150-150 150z m400 150h-120c-16 0-34 13-39 29l-31 93c-6 15-23 28-40 28h-340c-16 0-34-13-39-28l-31-94c-6-15-23-28-40-28h-120c-55 0-100-45-100-100v-450c0-55 45-100 100-100h800c55 0 100 45 100 100v450c0 55-45 100-100 100z m-400-550c-138 0-250 112-250 250 0 138 112 250 250 250 138 0 250-112 250-250 0-138-112-250-250-250z m365 380c-19 0-35 16-35 35 0 19 16 35 35 35 19 0 35-16 35-35 0-19-16-35-35-35z" transform="matrix(1 0 0 -1 0 850)"></path></svg>
            </button>
        </div>
        <div id="chart"></div>
        <script>
          // helper to trigger download of a Blob
          function triggerDownload(blob, filename) {{
            const url = URL.createObjectURL(blob);
            const a = document.createElement('a');
            a.href = url;
            a.download = filename;
            document.body.appendChild(a);
            a.click();
            document.body.removeChild(a);
            URL.revokeObjectURL(url);
          }}

          // SVG → file
          document.getElementById('download-svg').addEventListener('click', () => {{
            const svg = document.querySelector('#chart svg');
            const serializer = new XMLSerializer();
            let source = serializer.serializeToString(svg);
            // ensure namespace
            if (!source.match(/^<svg[^>]+xmlns="http:\/\/www\.w3\.org\/2000\/svg"/)) {{
              source = source.replace(/^<svg/, '<svg xmlns="http://www.w3.org/2000/svg"');
            }}
            const blob = new Blob([source], {{ type: 'image/svg+xml;charset=utf-8' }});
            triggerDownload(blob, 'chart.svg');
          }});

          // SVG → PNG
          document.getElementById('download-png').addEventListener('click', () => {{
            const svgEl = document.querySelector('#chart svg');
            // pull viewBox values
            const vb = svgEl.viewBox.baseVal;
            // options to capture full box (negative x/y)
            const options = {{
              left:   vb.x,
              top:    vb.y,
              width:  vb.width,
              height: vb.height,
              scale:  1,
              backgroundColor: null
            }};
            saveSvgAsPng(svgEl, 'chart.png', options);
          }});   
        </script>        
      <script>
        function renderNetwork(data, containerSelector = '#chart') {{
            // Specify the dimensions of the chart.
            const width = {width};
            const height = {height}-50;

            // Specify the color scale.
            const color = d3.scaleOrdinal(d3.schemeCategory10);

            // Clone links and nodes so re-running is deterministic.
            const links = data.links.map(d => ({{ ...d }}));
            const nodes = data.nodes.map(d => ({{ ...d }}));

            // Create a simulation with adjusted charge strength and slower decay for longer runtime.
            const simulation = d3
                .forceSimulation(nodes)
                .force(
                  "link",
                  d3
                    .forceLink(links)
                    .id((d) => d.id)
                    .strength(0.1)
                )
                .force("charge", d3.forceManyBody().strength(-50))
                .force("x", d3.forceX())
                .force("y", d3.forceY())
                .force(
                  "collision",
                  d3.forceCollide().radius((d) => (d.type === "gene" ? 5 : 10) + 5)
                )
                .alphaDecay(0.02);

            // Select the container and clear previous contents.
            const container = d3.select(containerSelector);
            container.selectAll('*').remove();

            // Create the SVG element within the container.
            const svg = container.append('svg')
                .attr('width', width)
                .attr('height', height)
                .attr('viewBox', [-width / 2, -height / 2, width, height])
                .attr('style', 'max-width: 100%; height: auto;');

            // Draw links
            const link = svg
                .append("g")
                .attr("stroke-opacity", 0.6)
                .selectAll("line")
                .data(links)
                .join("line")
                .attr("stroke", (d) => d.color)
                .attr("stroke-width", (d) => Math.sqrt(d.value))
                .attr("opacity", 1);

            // Draw nodes
            const node = svg
                .append("g")
                .selectAll("circle")
                .data(nodes)
                .join("circle")
                .attr("r", (d) => (d.type === "gene" ? 5 : 10))
                .attr("fill", (d) => d.color || "darkgray")
                .attr("stroke", "#fff")
                .attr("stroke-opacity", 0.6)
                .attr("stroke-width", 1)
                .attr("opacity", 1);

            // Labels, hidden by default
            const labels = svg
                .append("g")
                .attr("pointer-events", "none")
                .selectAll("text")
                .data(nodes)
                .join("text")
                .attr("text-anchor", "middle")
                .attr("font-size", 10)
                .attr("stroke", "#fff")
                .attr("stroke-opacity", 0.3)
                .attr("stroke-width", 0.3)
                .attr("font-family", "Arial")
                .attr("dy", (d) => (d.type === "gene" ? -5 : -10) - 2)
                .attr("opacity", 0)
                .text((d) => d.label);

            // Drag behavior
            node.call(
              d3.drag().on("start", dragstarted).on("drag", dragged).on("end", dragended)
            );

            // Hover interactions
            node
              .on("mouseover", function (event, d) {{
                // Pin hovered node
                d.fx = d.x;
                d.fy = d.y;

                // Highlight links and dim others
                link
                  .attr("stroke", (l) => l.color)
                  .attr("stroke-width", (l) =>
                    l.source.id === d.id || l.target.id === d.id ? 2 : Math.sqrt(l.value)
                  )
                  .attr("opacity", (l) =>
                    l.source.id === d.id || l.target.id === d.id ? 1 : 0.5
                  );

                // Identify connected nodes
                const connected = new Set();
                links.forEach((l) => {{
                  if (l.source.id === d.id) connected.add(l.target.id);
                  else if (l.target.id === d.id) connected.add(l.source.id);
                }});

                // Highlight nodes and dim others
                node
                  .attr("stroke", (n) =>
                    n.id === d.id || connected.has(n.id) ? n.color || "darkgray" : "#fff"
                  )
                  .attr("stroke-width", (n) =>
                    n.id === d.id || connected.has(n.id) ? 3 : 1
                  )
                  .attr("opacity", (n) =>
                    n.id === d.id || connected.has(n.id) ? 1 : 0.5
                  );

                // Show labels for highlighted nodes
                labels.attr("opacity", (n) =>
                  n.id === d.id || connected.has(n.id) ? 1 : 0
                );

                // Add temporary repulsion force
                simulation
                  .force(
                    "hoverRepel",
                    d3
                      .forceManyBody()
                      .strength((n) => (n.id === d.id ? -3000 : 0))
                      .distanceMax(300)
                  )
                  .alpha(0.5)
                  .restart();
              }})
              .on("mouseout", function (event, d) {{
                // Unpin
                d.fx = null;
                d.fy = null;

                // Reset links
                link
                  .attr("stroke", (l) => l.color)
                  .attr("stroke-width", (l) => Math.sqrt(l.value))
                  .attr("opacity", 1);

                // Reset nodes
                node.attr("stroke", "#fff").attr("stroke-width", 1).attr("opacity", 1);

                // Hide labels
                labels.attr("opacity", 0);

                // Remove repulsion force
                simulation.force("hoverRepel", null).alphaTarget(0);
              }});

            // Update positions on tick
            simulation.on("tick", () => {{
              link
                .attr("x1", (d) => d.source.x)
                .attr("y1", (d) => d.source.y)
                .attr("x2", (d) => d.target.x)
                .attr("y2", (d) => d.target.y);

              node.attr("cx", (d) => d.x).attr("cy", (d) => d.y);

              labels.attr("x", (d) => d.x).attr("y", (d) => d.y);
            }});

            function dragstarted(event) {{
              if (!event.active) simulation.alphaTarget(0.3).restart();
              event.subject.fx = event.subject.x;
              event.subject.fy = event.subject.y;
            }}
            function dragged(event) {{
              event.subject.fx = event.x;
              event.subject.fy = event.y;
            }}
            function dragended(event) {{
              if (!event.active) simulation.alphaTarget(0);
              event.subject.fx = null;
              event.subject.fy = null;
            }}

            // When this cell is re-run, stop the previous simulation. (This doesn’t
            // really matter since the target alpha is zero and the simulation will
            // stop naturally, but it’s a good practice.)

            return svg.node();
        }}
      </script>
    </body>
    <script>
      document.addEventListener('DOMContentLoaded', () => {{
        renderNetwork({data_json}, '#chart');
      }});
    </script>        
    </html>
    """
    st.components.v1.html(html, height=height, scrolling=True)
