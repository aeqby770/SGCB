"""Build revision documents; --documents-only reuses the formatted workbook."""
from pathlib import Path
import sys
import subprocess
import re
from docx import Document
from docx.shared import Inches, Pt, RGBColor
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from openpyxl import load_workbook
from openpyxl.styles import Font, PatternFill, Alignment

ROOT = Path(sys.argv[1])
MAN = ROOT / "manuscript"
RESPONSE = ROOT / "response"
TMP = ROOT / "tmp" / "documents"
TMP.mkdir(parents=True, exist_ok=True)
PANDOC = ROOT / ".tools" / "pandoc.exe"
OFFICE = ROOT / ".tools" / "LibreOffice" / "program" / "soffice.com"

default_reference = subprocess.run([str(PANDOC), "--print-default-data-file", "reference.docx"],
                                   check=True, capture_output=True).stdout
(TMP / "reference.docx").write_bytes(default_reference)
reference = Document(TMP / "reference.docx")
section = reference.sections[0]
section.page_width = Inches(8.2677)
section.page_height = Inches(11.6929)
section.top_margin = section.bottom_margin = Inches(.78)
section.left_margin = section.right_margin = Inches(.78)
styles = reference.styles
for name in ("Normal", "Body Text", "First Paragraph"):
    if name in styles:
        styles[name].font.name = "Times New Roman"
        styles[name].font.size = Pt(11)
        styles[name].paragraph_format.line_spacing = 1.22
        styles[name].paragraph_format.space_after = Pt(6)
for name, size in (("Heading 1",16),("Heading 2",13),("Heading 3",11.5),("Heading 4",11)):
    styles[name].font.name = "Arial"
    styles[name].font.size = Pt(size)
    styles[name].font.color.rgb = RGBColor(0,0,0)
    styles[name].paragraph_format.keep_with_next = True
    styles[name].paragraph_format.space_before = Pt(12)
    styles[name].paragraph_format.space_after = Pt(6)
styles["Caption"].font.name = "Times New Roman"
styles["Caption"].font.size = Pt(10)
reference.save(TMP / "reference.docx")

main_text = (MAN / "CSBJ_SGCB_manuscript.md").read_text(encoding="utf-8")
main_text += "\n\n" + (MAN / "tables.md").read_text(encoding="utf-8")
main_text += "\n\n" + (MAN / "figure_captions.md").read_text(encoding="utf-8")
(TMP / "main_assembled.md").write_text(main_text, encoding="utf-8")
jobs = [
    (TMP / "main_assembled.md", MAN / "CSBJ_SGCB_manuscript.docx", "SGCB · Revised manuscript"),
    (MAN / "CSBJ_SGCB_supplementary.md", MAN / "CSBJ_SGCB_supplementary.docx", "SGCB · Supplementary methods and results"),
    (RESPONSE / "Response_to_Editor_and_Reviewers.md", RESPONSE / "Response_to_Editor_and_Reviewers.docx", "SGCB · Response to the Editor and Reviewers")]
for source, target, label in jobs:
    subprocess.run([str(PANDOC), str(source), "--from=markdown+tex_math_dollars+pipe_tables", "--to=docx",
                    "--reference-doc=" + str(TMP / "reference.docx"), "--output=" + str(target)], check=True)
    doc = Document(target)
    for paragraph in doc.paragraphs:
        if paragraph.style.name.startswith("Heading"):
            ppr = paragraph._p.get_or_add_pPr()
            numpr = ppr.find(qn("w:numPr"))
            if numpr is None:
                numpr = OxmlElement("w:numPr")
                ppr.append(numpr)
            numid = numpr.find(qn("w:numId"))
            if numid is None:
                numid = OxmlElement("w:numId")
                numpr.append(numid)
            numid.set(qn("w:val"), "0")
        if paragraph.text in ("Tables", "Figure legends") or (
                source.name == "main_assembled.md" and paragraph.text.startswith("Table 3.")):
            paragraph.paragraph_format.page_break_before = True
        if paragraph.text.startswith("Figure "):
            paragraph.paragraph_format.keep_together = True
        if paragraph.text.startswith(("The primary qPCR effect threshold", "Brain contrasts each include")):
            paragraph.paragraph_format.line_spacing = 1.05
            paragraph.paragraph_format.space_after = Pt(3)
            for run in paragraph.runs:
                run.font.size = Pt(9)
    for section in doc.sections:
        header = section.header.paragraphs[0]
        header.text = label
        header.runs[0].font.name = "Arial"
        header.runs[0].font.size = Pt(8)
        header.runs[0].font.color.rgb = RGBColor.from_string("666666")
        footer = section.footer.paragraphs[0]
        footer.alignment = 2
        field = OxmlElement("w:fldSimple")
        field.set(qn("w:instr"), "PAGE")
        footer._p.append(field)
    for previous, paragraph in zip(doc.paragraphs, doc.paragraphs[1:]):
        if paragraph.text.startswith("Location."):
            previous.paragraph_format.keep_with_next = True
    for table in doc.tables:
        table.autofit = False
        cols = len(table.columns)
        widths = {3:[2.35,2.05,2.23],6:[1.0,1.85,.72,.72,1.0,1.0],5:[1.1,1.0,1.0,1.0,2.0],4:[1.0,1.8,1.8,2.03]}.get(cols)
        if table.cell(0, 0).text == "Aspect":
            widths = [1.05, 2.20, 3.38]
        if table.cell(0, 0).text == "Cohort":
            widths = [1.13,.8,.68,.5,1.45,2.07]
        if cols == 5 and table.cell(0, 1).text == "Target":
            widths = [1.2,1.7,1.0,1.35,1.38]
        if widths:
            for column, width in zip(table.columns, widths):
                column.width = Inches(width)
            for row in table.rows:
                for cell, width in zip(row.cells, widths):
                    cell.width = Inches(width)
        for i, row in enumerate(table.rows):
            for cell in row.cells:
                for p in cell.paragraphs:
                    p.paragraph_format.line_spacing = 1.05
                    p.paragraph_format.space_after = Pt(3)
                    for run in p.runs:
                        run.font.size = Pt(9)
                        run.font.name = "Times New Roman"
                        run.bold = (i == 0)
            if i == 0:
                repeat = OxmlElement("w:tblHeader")
                row._tr.get_or_add_trPr().append(repeat)
            no_split = OxmlElement("w:cantSplit")
            row._tr.get_or_add_trPr().append(no_split)
    if source.name == "main_assembled.md":
        for number in (1,2,3,4,5,6):
            doc.add_page_break()
            doc.add_heading(f"Figure {number}", level=2)
            doc.add_picture(str(MAN / "figures" / f"Figure_{number}.png"), width=Inches(6.65))
    if source.name == "CSBJ_SGCB_supplementary.md":
        for number in (1,2,3,4,5,6,7,8):
            doc.add_page_break()
            doc.add_heading(f"Figure S{number}", level=2)
            doc.add_picture(str(MAN / "figures" / f"Figure_S{number}.png"), width=Inches(6.65))
    if source.name == 'Response_to_Editor_and_Reviewers.md':
        doc.paragraphs[-1].paragraph_format.keep_together = True
    doc.core_properties.title = label
    doc.core_properties.author = "Hongjie Jiang; Anchi Che; Yanxiao Han; Yanshuo Han"
    doc.save(target)

if "--documents-only" not in sys.argv[2:]:
    book = load_workbook(MAN / "Supplementary_Tables.xlsx")
    for sheet in book.worksheets:
        sheet.freeze_panes = "A2"
        sheet.auto_filter.ref = sheet.dimensions
        sheet.sheet_view.showGridLines = False
        for cell in sheet[1]:
            cell.font = Font(name="Calibri", size=10, color="FFFFFF", bold=True)
            cell.fill = PatternFill("solid", fgColor="286080")
            cell.alignment = Alignment(wrap_text=True, vertical="center")
        sheet.row_dimensions[1].height = 32
        for column in sheet.columns:
            max_width = min(42, max(12, max(len(str(c.value or "")) for c in column) + 2))
            sheet.column_dimensions[column[0].column_letter].width = max_width
            for cell in column[1:]:
                cell.font = Font(name="Calibri", size=10)
                if isinstance(cell.value, float):
                    header_name = str(column[0].value)
                    cell.number_format = "0.000E+00" if re.search(r"(^p$|^q$|p[._]?value|q[._]?value|_[pq]$)", header_name, re.I) else "0.0000"
        sheet.print_title_rows = "1:1"
        sheet.sheet_properties.pageSetUpPr.fitToPage = True
        sheet.page_setup.orientation = "landscape"
        sheet.page_setup.paperSize = sheet.PAPERSIZE_A3
        sheet.page_setup.fitToWidth = 1
        sheet.page_setup.fitToHeight = 0
    book.save(MAN / "Supplementary_Tables.xlsx")

profile = (TMP / "libreoffice_profile").as_uri()
for _, target, _ in jobs:
    completed = subprocess.run([str(OFFICE), "-env:UserInstallation=" + profile, "--headless",
        "--convert-to", "pdf", "--outdir", str(target.parent), str(target)],
        check=True, capture_output=True, text=True, encoding="utf-8", errors="replace",
        creationflags=subprocess.CREATE_NO_WINDOW)
    print("Rendered", target.with_suffix(".pdf").name)
print("Built three editable DOCX documents and their PDFs.")
