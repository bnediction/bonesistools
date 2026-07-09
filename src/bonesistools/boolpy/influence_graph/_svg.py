#!/usr/bin/env python

from __future__ import annotations

import re
from typing import Optional, Union

SvgLength = Union[str, int, float]


def scale_svg(
    svg: str,
    width: Optional[SvgLength] = None,
    height: Optional[SvgLength] = None,
) -> str:
    """
    Return SVG text with optional root width and height attributes.
    """

    if width is None and height is None:
        return svg

    start = svg.find("<svg")

    if start == -1:
        return svg

    end = svg.find(">", start)

    if end == -1:
        return svg

    root = svg[start:end]

    for attribute in ("width", "height"):
        root = re.sub(
            rf"\s{attribute}\s*=\s*(['\"]).*?\1",
            "",
            root,
            count=1,
        )

    attributes = []

    if width is not None:
        attributes.append(f'width="{width}"')

    if height is not None:
        attributes.append(f'height="{height}"')

    if attributes:
        root = f"{root} {' '.join(attributes)}"

    return f"{svg[:start]}{root}{svg[end:]}"
