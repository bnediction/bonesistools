#!/usr/bin/env python

from typing import Sequence

from itertools import cycle

def rgb(color: list):
    return list(map(lambda x: x/255, color))

def rgb2hex(rgb: Sequence):
    r, g, b = rgb
    if isinstance(r, float) and isinstance(g, float) and isinstance(b, float):
        r = round(r*255) if 0<r<1 else int(r)
        g = round(g*255) if 0<g<1 else int(g)
        b = round(b*255) if 0<b<1 else int(b)
    return "#{:02x}{:02x}{:02x}".format(r, g, b)

black       = rgb([  0,   0,   0])
white       = rgb([255, 255, 255])
blue        = rgb([  0,  20, 255])
red         = rgb([255,  80,  50])
green       = rgb([ 20, 200,  80])
violet      = rgb([255,  51, 255])
lightgreen  = rgb([ 20, 250,  80])
coral       = rgb([255, 127,  80])
yellow      = rgb([255, 255,   0])
darkyellow  = rgb([204, 204,   0])
lightyellow = rgb([128, 128,   0])
darkorange  = rgb([255, 140,   0])
lightorange = rgb([255, 165,  90])
limegreen   = rgb([ 50, 255,  50])
pink        = rgb([255, 182, 193])
orchid      = rgb([218, 112, 214])
magenta     = rgb([255,   0, 255])
purple      = rgb([128,   0, 128])
indigo      = rgb([ 75,   0, 130])
slateblue   = rgb([ 71,  60, 139])
lightgray   = rgb([211, 211, 211])
gray        = rgb([112, 128, 144])
darkgreen   = rgb([  0, 100,   0])
gold        = rgb([238, 201,   0])
orange      = rgb([255, 165,   0])
salmon      = rgb([198, 113, 113])
maroon      = rgb([128,   0,   0])
beet        = rgb([142,  56, 142])
teal        = rgb([ 56, 142, 142])
olive       = rgb([142, 142,  56])
navy        = rgb([  0,   0, 128])
darkblue    = rgb([  0,   0, 139])
skyblue     = rgb([135, 206, 235])

COLORS = [
    blue,
    red,
    green,
    orange,
    purple,
    skyblue,
    teal,
    pink,
    violet,
    darkblue,
    magenta,
    darkgreen,
    darkorange,
    gray,
    maroon,
    olive,
    orchid,
    beet,
    indigo,
    gold,
    navy,
    salmon
]

LIGHT_COLORS = [
    skyblue,
    red,
    green,
    orange,
    orchid,
    gray,
    teal,
    magenta,
    violet,
    olive,
    beet,
    indigo,
    gold,
    navy,
    salmon
]

color_cycle = cycle(COLORS)