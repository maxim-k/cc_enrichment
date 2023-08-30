function onKeyUp(event) {
  Streamlit.setComponentValue(event.target.value)
}

const delay = (callback, wait) => {
  let timeoutId = null;
  return (...args) => {
    window.clearTimeout(timeoutId);
    timeoutId = window.setTimeout(() => {
      callback.apply(null, args);
    }, wait);
  };
}

/**
 * The component's render function. This will be called immediately after
 * the component is initially loaded, and then again every time the
 * component gets new data from Python.
 */
function onRender(event) {
  // Get the RenderData from the event
  if (!window.rendered) {
    const {
      label,
      value,
      delay: delay_time,
      max_chars,
      height,
      type,
      placeholder,
      disabled,
      label_visibility
    } = event.detail.args;

    const input = document.getElementById("text_input");
    const label_el = document.getElementById("label")
    const root = document.getElementById("root")
    let frame_height;

    if (label_el) {
      label_el.innerText = label
    }

    if (value && !input.value) {
      input.value = value
    }

    input.type = "text"

    if (max_chars) {
      input.maxLength = max_chars
    }

    if (height === undefined || height === null) {
      frame_height = 129;
    } else {
      // If 'height' is set, set FrameHeight to the bigger number of the value or 73 (height of a one line input)
      frame_height = Math.max(height, 73);
    }

    if (placeholder) {
      input.placeholder = placeholder
    }

    if (disabled) {
      input.disabled = true
      label.disabled = true
      root.classList.add("disabled")
    }

    if (label_visibility == "hidden") {
      root.classList.add("label-hidden")
    }
    else if (label_visibility == "collapsed") {
      root.classList.add("label-collapsed")
      frame_height = frame_height - 29
    }

    if (delay_time > 0) { // is false if delay_time is 0 or undefined
      input.onkeyup = delay(onKeyUp, delay_time)
    }
    else {
      input.onkeyup = onKeyUp
    }

    Streamlit.setFrameHeight(frame_height)
    window.rendered = true
  }
}

Streamlit.events.addEventListener(Streamlit.RENDER_EVENT, onRender)
Streamlit.setComponentReady()