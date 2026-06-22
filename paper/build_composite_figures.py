
from __future__ import annotations

from pathlib import Path

from PIL import Image, ImageChops, ImageDraw, ImageFont

ROOT = Path(__file__).resolve().parent
IMG = ROOT / "images"
WHITE = (255, 255, 255)
INK = (20, 20, 20)

def _font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont | ImageFont.ImageFont:
    names = ["arialbd.ttf" if bold else "arial.ttf", "Arial Bold.ttf" if bold else "Arial.ttf"]
    for name in names:
        for base in [Path("C:/Windows/Fonts"), Path("/usr/share/fonts/truetype/dejavu")]:
            path = base / name
            if path.exists():
                return ImageFont.truetype(str(path), size)
    return ImageFont.load_default()

LABEL_FONT = _font(30, bold=True)

def trim_white(im: Image.Image, pad: int = 18, threshold: int = 248) -> Image.Image:

    rgb = im.convert("RGB")
    mask = Image.new("RGB", rgb.size, WHITE)
    diff = ImageChops.difference(rgb, mask).convert("L")
    diff = diff.point(lambda p: 255 if p > 255 - threshold else 0)
    box = diff.getbbox()
    if not box:
        return rgb
    left, top, right, bottom = box
    left = max(0, left - pad)
    top = max(0, top - pad)
    right = min(rgb.width, right + pad)
    bottom = min(rgb.height, bottom + pad)
    return rgb.crop((left, top, right, bottom))

def load(name: str) -> Image.Image:
    return trim_white(Image.open(IMG / name))

def crop_fraction(im: Image.Image, box: tuple[float, float, float, float]) -> Image.Image:
    left, top, right, bottom = box
    return im.crop(
        (
            int(left * im.width),
            int(top * im.height),
            int(right * im.width),
            int(bottom * im.height),
        )
    )

def fit(im: Image.Image, max_w: int, max_h: int) -> Image.Image:
    out = im.copy()
    out.thumbnail((max_w, max_h), Image.Resampling.LANCZOS)
    return out

def label_panel(canvas: Image.Image, xy: tuple[int, int], label: str) -> None:
    draw = ImageDraw.Draw(canvas)
    x, y = xy
    bbox = draw.textbbox((0, 0), label, font=LABEL_FONT)
    pad_x, pad_y = 12, 7
    rect = (x, y, x + bbox[2] + 2 * pad_x, y + bbox[3] + 2 * pad_y)
    draw.rounded_rectangle(rect, radius=8, fill=(255, 255, 255), outline=(210, 210, 210))
    draw.text((x + pad_x, y + pad_y), label, fill=INK, font=LABEL_FONT)

def stack_vertical(
    panels: list[tuple[str, str]],
    output: str,
    width: int = 1800,
    max_panel_h: int = 720,
    gap: int = 36,
    margin: int = 45,
) -> None:
    rendered = [(fit(load(name), width - 2 * margin, max_panel_h), label) for name, label in panels]
    height = 2 * margin + sum(im.height for im, _ in rendered) + gap * (len(rendered) - 1)
    canvas = Image.new("RGB", (width, height), WHITE)
    y = margin
    for im, label in rendered:
        x = (width - im.width) // 2
        canvas.paste(im, (x, y))
        label_panel(canvas, (x + 14, y + 14), label)
        y += im.height + gap
    canvas.save(IMG / output, optimize=True)

def row_equal_height(
    panels: list[tuple[str, str]],
    output: str,
    height: int = 560,
    gap: int = 34,
    margin: int = 45,
) -> None:
    rendered = []
    for name, label in panels:
        im = load(name)
        scale = height / im.height
        rendered.append((im.resize((int(im.width * scale), height), Image.Resampling.LANCZOS), label))
    width = 2 * margin + sum(im.width for im, _ in rendered) + gap * (len(rendered) - 1)
    canvas = Image.new("RGB", (width, height + 2 * margin), WHITE)
    x = margin
    for im, label in rendered:
        y = margin
        canvas.paste(im, (x, y))
        label_panel(canvas, (x + 12, y + 12), label)
        x += im.width + gap
    canvas.save(IMG / output, optimize=True)

def grid(
    panels: list[tuple[str, str]],
    output: str,
    cols: int,
    panel_w: int = 720,
    panel_h: int = 500,
    gap: int = 30,
    margin: int = 38,
) -> None:
    rows = (len(panels) + cols - 1) // cols
    width = 2 * margin + cols * panel_w + (cols - 1) * gap
    height = 2 * margin + rows * panel_h + (rows - 1) * gap
    canvas = Image.new("RGB", (width, height), WHITE)
    for row in range(rows):
        row_panels = panels[row * cols : (row + 1) * cols]
        row_width = len(row_panels) * panel_w + max(0, len(row_panels) - 1) * gap
        x0 = (width - row_width) // 2
        for col, (name, label) in enumerate(row_panels):
            x = x0 + col * (panel_w + gap)
            y = margin + row * (panel_h + gap)
            im = fit(load(name), panel_w, panel_h)
            canvas.paste(im, (x + (panel_w - im.width) // 2, y + (panel_h - im.height) // 2))
            label_panel(canvas, (x + 10, y + 10), label)
    canvas.save(IMG / output, optimize=True)

def nooru_mesh_layout(output: str) -> None:
    width = 1650
    margin = 34
    gap = 42
    panel_w = (width - 2 * margin - gap) // 2
    panel_h = 650

    bc = fit(load("nooru_BC_2D.png"), panel_w, 430)
    mesh = fit(load("nooru_mesh_3D.png"), panel_w, panel_h)

    crack = fit(crop_fraction(load("Exp_noor.png"), (0.0, 0.12, 1.0, 0.76)), panel_w, 170)

    height = 2 * margin + panel_h
    canvas = Image.new("RGB", (width, height), WHITE)

    x = margin
    y = margin
    canvas.paste(bc, (x + (panel_w - bc.width) // 2, y))
    canvas.paste(crack, (x + (panel_w - crack.width) // 2, y + panel_h - crack.height))
    label_panel(canvas, (x + 10, y + 10), "(a)")

    x = margin + panel_w + gap
    canvas.paste(mesh, (x + (panel_w - mesh.width) // 2, y + (panel_h - mesh.height) // 2))
    label_panel(canvas, (x + 10, y + 10), "(b)")
    canvas.save(IMG / output, optimize=True)

def three_pb_layout(output: str) -> None:

    grid(
        [
            ("fig_mesh.png", "(a)"),
            ("abaqus_fig_damage_last_step.png", "(b)"),
            ("fig_damage_peak.png", "(c)"),
            ("fig_damage_postpeak.png", "(d)"),
        ],
        output,
        cols=2,
        panel_w=790,
        panel_h=405,
        gap=22,
        margin=26,
    )

def main() -> None:
    stack_vertical(
        [
            ("fig_mesh.png", "(a)"),
            ("abaqus_fig_damage_last_step.png", "(b)"),
        ],
        "fig_b1_mesh_abq.png",
        width=1650,
        max_panel_h=520,
        gap=24,
        margin=30,
    )
    stack_vertical(
        [
            ("fig_damage_peak.png", "(a)"),
            ("fig_damage_postpeak.png", "(b)"),
        ],
        "fig_b1_damage.png",
        width=1650,
        max_panel_h=390,
        gap=22,
        margin=28,
    )
    row_equal_height(
        [
            ("load_cmod_comparison.png", "(a)"),
            ("time_comparison_bar.png", "(b)"),
        ],
        "fig_b1_results.png",
        height=430,
        gap=28,
        margin=30,
    )
    three_pb_layout("fig_b1_compact.png")
    nooru_mesh_layout("fig_b2_mesh.png")
    grid(
        [
            ("nooru_inc_0029.png", "(a)"),
            ("nooru_inc_0037.png", "(b)"),
            ("nooru_inc_0041.png", "(c)"),
            ("nooru_inc_0053.png", "(d)"),
            ("nooru_inc_0075.png", "(e)"),
            ("nooru_inc_0127.png", "(f)"),
            ("nooru_inc_0275.png", "(g)"),
            ("nooru_inc_0407.png", "(h)"),
            ("nooru_inc_0900.png", "(i)"),
        ],
        "nooru_damage_evolution_3x3.png",
        cols=3,
        panel_w=510,
        panel_h=330,
        gap=18,
        margin=20,
    )
    row_equal_height(
        [
            ("torsion.png", "(a)"),
            ("Exp_torsion.png", "(b)"),
        ],
        "fig_b3_mesh.png",
        height=430,
        gap=32,
        margin=32,
    )
    grid(
        [
            ("Job-1_StaticFast_mod_vm_LIVE_snap_inc_0001_theta_2_143e-05.png", "(a)"),
            ("Job-1_StaticFast_mod_vm_LIVE_snap_inc_0021_theta_4_500e-04.png", "(b)"),
            ("Job-1_StaticFast_mod_vm_LIVE_snap_inc_0041_theta_8_786e-04.png", "(c)"),
            ("Job-1_StaticFast_mod_vm_LIVE_snap_inc_0061_theta_1_307e-03.png", "(d)"),
            ("Job-1_StaticFast_mod_vm_LIVE_snap_inc_0080_theta_1_714e-03.png", "(e)"),
            ("Job-1_StaticFast_mod_vm_LIVE_snap_inc_0100_theta_2_143e-03.png", "(f)"),
            ("Job-1_StaticFast_mod_vm_LIVE_snap_inc_0120_theta_2_571e-03.png", "(g)"),
            ("Job-1_StaticFast_mod_vm_LIVE_snap_inc_0140_theta_3_000e-03.png", "(h)"),
        ],
        "fig_b3_damage_evolution.png",
        cols=3,
        panel_w=520,
        panel_h=300,
        gap=18,
        margin=20,
    )

if __name__ == "__main__":
    main()
