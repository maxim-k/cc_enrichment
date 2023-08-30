from pathlib import Path
from typing import Any, Callable, Dict, Optional, Tuple

import streamlit as st
import streamlit.components.v1 as components

frontend_dir = (Path(__file__).parent / "frontend").absolute()
_component_func = components.declare_component(
    "textarea_onchange", path=str(frontend_dir)
)


def textarea_onchange(
    label: str,
    value: str = "",
    max_chars: Optional[int] = None,
    height: Optional[int] = None,
    key: Optional[str] = None,
    type: str = "default",
    delay: Optional[int] = None,
    on_change: Optional[Callable] = None,
    args: Optional[Tuple[Any, ...]] = None,
    kwargs: Optional[Dict[str, Any]] = None,
    *,
    placeholder: str = "",
    disabled: bool = False,
    label_visibility: str = "visible",
):
    """
    This component is based on https://github.com/blackary/streamlit-keyup
    Instead of <input/> and st.text_input it uses <textarea/> and st.text_area,
    as in Stremlit they are nearly identical, and you can replace
    st.text_input with st.text_area with height of 73.

    Generate a text input that renders on keyup.

    Delay means that it will wait at least the specified amount of milliseconds
    before updating the value. This is useful for preventing excessive updates
    when the user is typing. Since the input updating will cause the app to rerun,
    if you are having performance issues, you should consider setting a delay
    value.

    on_change is a callback function that will be called when the value changes.

    args and kwargs are optional arguments which are passed to the on_change callback
    function
    """

    if key is None:
        key = "textarea_onchange_" + label

    component_value = _component_func(
        label=label,
        value=value,
        key=key,
        default=value,
        max_chars=max_chars,
        height=height,
        type=type,
        delay=delay,
        placeholder=placeholder,
        disabled=disabled,
        label_visibility=label_visibility,
    )

    if on_change is not None:
        if "__previous_values__" not in st.session_state:
            st.session_state["__previous_values__"] = {}

        if component_value != st.session_state["__previous_values__"].get(key, value):
            st.session_state["__previous_values__"][key] = component_value

            if on_change:
                if args is None:
                    args = ()
                if kwargs is None:
                    kwargs = {}
                on_change(*args, **kwargs)

    return component_value


def main():
    return


if __name__ == "__main__":
    main()
