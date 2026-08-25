"""Shared builders for the dice scenes.

The entry-point scripts (dice-light-area.py, dice-light-spot.py,
dice-sentence.py) import this module and only decide composition: which
dice, which light, which camera. Run them from the repo root; Python puts
scenes/ on the import path automatically when the script lives there.

A die is a CSG chain: intersection(cube, sphere r=1.4) trims the corners,
then 21 pip spheres (hemispherical dimples, white material) are subtracted
one difference at a time. Materials sit on the CSG leaves — a material on a
CSG node would recolor the whole subtree, pips included.
"""

import json
import math
import random
import sys

PIP_R = 0.2
# Corner-trimming sphere. Edges sit at sqrt(2)~1.414 from the center, corners
# at sqrt(3)~1.732, so r=1.4 rounds the corners and just shaves the edges.
TRIM_R = 1.4

# Pip layouts per face value: (u, v) offsets in the face plane.
O = 0.5
LAYOUT = {
    1: [(0, 0)],
    2: [(-O, -O), (O, O)],
    3: [(-O, -O), (0, 0), (O, O)],
    4: [(-O, -O), (-O, O), (O, -O), (O, O)],
    5: [(-O, -O), (-O, O), (0, 0), (O, -O), (O, O)],
    6: [(-O, -O), (-O, 0), (-O, O), (O, -O), (O, 0), (O, O)],
}

# face -> (value, axis, sign); opposite faces sum to 7.
FACES = [
    (1, "y", +1), (6, "y", -1),
    (3, "x", +1), (4, "x", -1),
    (2, "z", +1), (5, "z", -1),
]

PIP_MATERIAL = {
    "color": [0.96, 0.94, 0.88],
    "diffuse": 0.8, "ambient": 0.15, "specular": 0.4, "shininess": 100,
}

# Vibrant (light, dark) marble hue pairs, bowling-ball style.
MAGENTA = ([0.93, 0.08, 0.45], [0.45, 0.02, 0.28])
TEAL = ([0.0, 0.78, 0.82], [0.02, 0.32, 0.45])
ORANGE = ([1.0, 0.55, 0.05], [0.75, 0.18, 0.05])
LIME = ([0.55, 0.92, 0.1], [0.22, 0.5, 0.05])
VIOLET = ([0.6, 0.25, 0.95], [0.3, 0.08, 0.5])
GOLD = ([1.0, 0.8, 0.12], [0.62, 0.42, 0.05])
PALETTE = [MAGENTA, TEAL, ORANGE, LIME, VIOLET, GOLD]


def checkers(a):
    """The room's checkers pattern; `a` sets each surface's light-square tint."""
    return {
        "type": "checkers",
        "a": a, "b": [0.14, 0.13, 0.16],
        "transform": [{"scale": [2, 2, 2]}],
    }


def marble(a, b):
    """Bowling-ball swirl: perturbed stripes in one vibrant hue."""
    return {
        "type": "perturbed",
        "amplitude": 1.1,
        "frequency": 0.9,
        "seed": 11,
        "pattern": {
            "type": "stripe", "a": a, "b": b,
            "transform": [{"scale": [0.35, 0.35, 0.35]}, {"rotate_z": 0.9}],
        },
    }


def body_material(a, b):
    return {
        "pattern": marble(a, b),
        "diffuse": 0.75, "ambient": 0.2,
        "specular": 1.0, "shininess": 300, "reflective": 0.08,
    }


def pip(axis, sign, u, v):
    """A dimple sphere centered on the face at (u, v)."""
    center = {
        "x": [sign * 1.0, u, v],
        "y": [u, sign * 1.0, v],
        "z": [u, v, sign * 1.0],
    }[axis]
    return {
        "type": "sphere",
        "material": dict(PIP_MATERIAL),
        "transform": [{"scale": [PIP_R, PIP_R, PIP_R]}, {"translate": center}],
    }


def die(colors, transform):
    """One die (unit half-width) as a CSG tree, with `transform` applied last."""
    a, b = colors
    node = {
        "type": "csg", "operation": "intersection",
        "left": {"type": "cube", "material": body_material(a, b)},
        "right": {
            "type": "sphere", "material": body_material(a, b),
            "transform": [{"scale": [TRIM_R, TRIM_R, TRIM_R]}],
        },
    }
    for value, axis, sign in FACES:
        for u, v in LAYOUT[value]:
            node = {
                "type": "csg", "operation": "difference",
                "left": node,
                "right": pip(axis, sign, u, v),
            }
    node["transform"] = transform
    return node


def room(back_z=7, right_x=7):
    """Two checkered walls meeting in a corner, and a glossy checkered floor."""
    return [
        {
            # back wall, checkered in a pale mint tint
            "type": "plane",
            "material": {
                "pattern": checkers([0.85, 0.93, 0.88]),
                "diffuse": 0.8, "ambient": 0.15, "specular": 0.05,
            },
            "transform": [{"rotate_x": math.pi / 2},
                          {"translate": [0, 0, back_z]}],
        },
        {
            # right wall, a deeper mint so the corner reads
            "type": "plane",
            "material": {
                "pattern": checkers([0.75, 0.93, 0.88]),
                "diffuse": 0.8, "ambient": 0.15, "specular": 0.05,
            },
            "transform": [{"rotate_z": math.pi / 2},
                          {"translate": [right_x, 0, 0]}],
        },
        {
            # floor: warm cream checkers, glossy and slightly reflective
            "type": "plane",
            "material": {
                "pattern": checkers([0.95, 0.93, 0.88]),
                "diffuse": 0.8, "ambient": 0.2,
                "specular": 0.3, "reflective": 0.12,
            },
        },
    ]


def camera(width, height, from_, to, field_of_view=0.8, antialias=3):
    cam = {
        "width": width, "height": height, "field_of_view": field_of_view,
        "from": from_, "to": to, "up": [0, 1, 0],
    }
    if antialias:
        cam["antialias"] = {"n": antialias}
    return cam


def area_light():
    """Rectangular area light up and to the camera's left: soft shadows."""
    return {
        "type": "area",
        "corner": [-8, 12, -10], "uvec": [4, 0, 0], "usteps": 4,
        "vvec": [0, 0, 4], "vsteps": 4, "intensity": [1, 1, 1],
    }


def spot_light(position=(-3.0, 10.0, -4.0), target=(0.1, 2.0, 0.3),
               inner_angle=0.22, outer_angle=0.42, intensity=1.2):
    """Spotlight aimed at `target`: full intensity inside inner_angle, fading
    to nothing at outer_angle, so the subject sits in a pool of light."""
    return {
        "type": "spot",
        "position": list(position), "target": list(target),
        "intensity": [intensity] * 3,
        "inner_angle": inner_angle, "outer_angle": outer_angle,
    }


# 3x5 dot-matrix font: 5 rows (top first) of 3 cells; '#' becomes a die.
FONT = {
    "A": [".#.", "#.#", "###", "#.#", "#.#"],
    "B": ["##.", "#.#", "##.", "#.#", "##."],
    "C": [".##", "#..", "#..", "#..", ".##"],
    "D": ["##.", "#.#", "#.#", "#.#", "##."],
    "E": ["###", "#..", "##.", "#..", "###"],
    "F": ["###", "#..", "##.", "#..", "#.."],
    "G": [".##", "#..", "#.#", "#.#", ".##"],
    "H": ["#.#", "#.#", "###", "#.#", "#.#"],
    "I": ["###", ".#.", ".#.", ".#.", "###"],
    "J": ["..#", "..#", "..#", "#.#", ".#."],
    "K": ["#.#", "#.#", "##.", "#.#", "#.#"],
    "L": ["#..", "#..", "#..", "#..", "###"],
    "M": ["#.#", "###", "###", "#.#", "#.#"],
    "N": ["##.", "#.#", "#.#", "#.#", "#.#"],
    "O": ["###", "#.#", "#.#", "#.#", "###"],
    "P": ["###", "#.#", "###", "#..", "#.."],
    "Q": ["###", "#.#", "#.#", "###", "..#"],
    "R": ["###", "#.#", "##.", "#.#", "#.#"],
    "S": [".##", "#..", ".#.", "..#", "##."],
    "T": ["###", ".#.", ".#.", ".#.", ".#."],
    "U": ["#.#", "#.#", "#.#", "#.#", "###"],
    "V": ["#.#", "#.#", "#.#", "#.#", ".#."],
    "W": ["#.#", "#.#", "###", "###", "#.#"],
    "X": ["#.#", "#.#", ".#.", "#.#", "#.#"],
    "Y": ["#.#", "#.#", ".#.", ".#.", ".#."],
    "Z": ["###", "..#", ".#.", "#..", "###"],
    "0": ["###", "#.#", "#.#", "#.#", "###"],
    "1": [".#.", "##.", ".#.", ".#.", "###"],
    "2": ["##.", "..#", ".#.", "#..", "###"],
    "3": ["###", "..#", ".#.", "..#", "###"],
    "4": ["#.#", "#.#", "###", "..#", "..#"],
    "5": ["###", "#..", "##.", "..#", "##."],
    "6": [".##", "#..", "###", "#.#", "###"],
    "7": ["###", "..#", ".#.", ".#.", ".#."],
    "8": ["###", "#.#", "###", "#.#", "###"],
    "9": ["###", "#.#", "###", "..#", "##."],
    "!": [".#.", ".#.", ".#.", "...", ".#."],
    "-": ["...", "...", "###", "...", "..."],
    ".": ["...", "...", "...", "...", ".#."],
    " ": ["...", "...", "...", "...", "..."],
}


def letter_dice(ch, x_left, z, colors, scale, spacing, rng):
    """Dice for one glyph; x_left is the left edge, letters stand on y=0."""
    rows = FONT[ch]
    dice = []
    for r, row in enumerate(rows):
        for c, cell in enumerate(row):
            if cell != "#":
                continue
            x = x_left + (c + 0.5) * spacing
            y = scale + (len(rows) - 1 - r) * spacing
            dice.append(die(colors, [
                # quarter-turns show a random face; small jitter reads as
                # hand-placed rather than machined
                {"rotate_x": rng.randrange(4) * math.pi / 2},
                {"rotate_y": rng.randrange(4) * math.pi / 2},
                {"rotate_z": rng.uniform(-0.05, 0.05)},
                {"rotate_y": rng.uniform(-0.06, 0.06)},
                {"scale": [scale, scale, scale]},
                {"translate": [x, y, z]},
            ]))
    return dice


def sentence_dice(text, z=5.0, scale=0.35, gap=0.1, palette=None, seed=5):
    """All dice for `text` (one group per letter, colors cycling through the
    palette), centered on x=0 and standing on the floor at depth `z`.
    Returns (shapes, total_width). Raises KeyError for unsupported glyphs."""
    palette = palette or PALETTE
    text = text.upper()
    unsupported = sorted({c for c in text if c not in FONT})
    if unsupported:
        raise KeyError(
            f"unsupported: {unsupported}; glyphs: {''.join(sorted(FONT))}"
        )
    rng = random.Random(seed)  # fixed seed: regenerating gives the same scene
    spacing = 2 * scale + gap  # center-to-center distance of adjacent dice
    letter_w = 3 * spacing
    advance = letter_w + spacing  # one empty column between letters
    total_w = len(text) * advance - spacing
    shapes = []
    color_i = 0
    for i, ch in enumerate(text):
        x_left = -total_w / 2 + i * advance
        if ch == " ":
            continue
        colors = palette[color_i % len(palette)]
        color_i += 1
        letter = letter_dice(ch, x_left, z, colors, scale, spacing, rng)
        if letter:
            # one group per letter: root-level groups are what the loader's
            # bvh_threshold subdivides
            shapes.append({"type": "group", "children": letter})
    return shapes, total_w


def write_scene(scene, out):
    with open(out, "w") as f:
        json.dump(scene, f, indent=2)
        f.write("\n")
    width = scene["camera"]["width"]
    height = scene["camera"]["height"]
    print(f"wrote {out} ({width}x{height})")


def main(out, scene_fn, default_size=(3456, 2234)):
    """Entry point for the fixed-composition scripts: parse an optional
    `width height` from argv, build via scene_fn(width, height), write."""
    width, height = default_size
    if len(sys.argv) == 3:
        width, height = int(sys.argv[1]), int(sys.argv[2])
    write_scene(scene_fn(width, height), out)
