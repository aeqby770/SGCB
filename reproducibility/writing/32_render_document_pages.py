"""Render every document page and contact sheets for visual inspection."""
from pathlib import Path
import sys
import fitz
from PIL import Image, ImageDraw

root = Path(sys.argv[1])
output = root / "tmp" / "documents" / "qa"
output.mkdir(parents=True, exist_ok=True)
files = [root / "manuscript" / "CSBJ_SGCB_manuscript.pdf",
         root / "manuscript" / "CSBJ_SGCB_supplementary.pdf",
         root / "response" / "Response_to_Editor_and_Reviewers.pdf",
         root / "manuscript" / "CSBJ_SGCB_manuscript_marked.pdf",
         root / "cover_letter" / "CSBJ_cover_letter.pdf"]
for path in files:
    pdf = fitz.open(path)
    print(path.name, "pages=", len(pdf), "text_characters=", sum(len(p.get_text()) for p in pdf))
    for start in range(0, len(pdf), 12):
        canvas = Image.new("RGB", (1440, 1560), "#dfe3e6")
        draw = ImageDraw.Draw(canvas)
        for offset, page in enumerate(list(pdf)[start:start+12]):
            pix = page.get_pixmap(matrix=fitz.Matrix(0.55, 0.55), alpha=False)
            tile = Image.frombytes("RGB", [pix.width,pix.height], pix.samples)
            tile.thumbnail((344, 483))
            x = (offset % 4) * 360 + 8
            y = (offset // 4) * 520 + 23
            canvas.paste(tile, (x,y))
            draw.text((x,y-17),f"Page {start+offset+1}",fill="black")
        canvas.save(output / f"{path.stem}_contact_{start//12+1}.png")
    for i,page in enumerate(pdf):
        page.get_pixmap(matrix=fitz.Matrix(1.4,1.4),alpha=False).save(output / f"{path.stem}_page_{i+1}.png")
        print("detail",path.stem,i+1)
