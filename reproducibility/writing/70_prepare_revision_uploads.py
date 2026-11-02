"""Prepare journal documents and separate publication and reproduction archives."""
from pathlib import Path
import sys
import subprocess
import os
from datetime import date
from docx import Document

ROOT = Path(sys.argv[1]).resolve()
ORIGINAL = ROOT / "original_submission"
COVER = ROOT / "cover_letter"
MAN = ROOT / "manuscript"
VERSION = (ROOT / "software" / "SGCB" / "DESCRIPTION").read_text(encoding="utf-8").split("Version: ",1)[1].splitlines()[0]
TMP = ROOT / "tmp" / "documents"
OFFICE = ROOT / ".tools" / "LibreOffice" / "program" / "soffice.com"
COVER.mkdir(exist_ok=True)
TMP.mkdir(parents=True, exist_ok=True)

letter = Document(ORIGINAL / "cover_letter" / "CSBJ_cover_letter.docx")
replacements = {
    0: date.today().strftime("%d %B %Y").lstrip("0"),
    2: "Re: Major revision of CSBJ-D-26-00876",
    4: 'We resubmit our Research Article under CSBJ-D-26-00876 following the decision of 8 October 2026. The submitted version, SGCB 0.4.1, was entitled "SGCB: Distribution-Aware Differential Analysis Across Omics Modalities." The revised title is "SGCB: zero-aware analysis of mean and relative variability in omics data."',
    5: "In response to the statistical comments, we replaced the submitted log-abundance moderated-t test with unconditional arithmetic-mean inference and replaced the gene-wise BF/Cauchy DV calibration with paired zero-aware population inference. We corrected the error-rate definitions and withdrew the original GTEx brain leading-yield claim. Known-truth controls and matched Gamma/lognormal restrictions now quantify the contribution and cost of GG modeling.",
    6: f"We also corrected discovery/follow-up separation in KIRC, recomputed the gene sets and patient scores, added external assessment, and developed purified-monocyte SLE and metabolic-liver cases with available covariate sensitivities. The response opens with a 0.4.1-to-{VERSION} change summary, then reproduces all reviewer comments with the corresponding action, evidence and location, including reasons for scope choices. Clean and marked manuscripts, six main and eight supplementary figures, complete tables and SGCB {VERSION} installation/reproduction materials accompany the revision.",
    7: "Authorship, affiliations, funding and competing-interest statements are unchanged. The work is original and is not under consideration elsewhere. During the revision of this manuscript, the authors used OpenAI Codex to assist software development, data analysis, and manuscript preparation. The authors reviewed and edited the output and take full responsibility for the content.",
}
for index, value in replacements.items():
    letter.paragraphs[index].text = value
letter.paragraphs[2].runs[0].bold = True
letter.paragraphs[8].insert_paragraph_before(
    "Author contacts: Hongjie Jiang, jiangxiantian2@gmail.com; Anchi Che, cheanchi0712@163.com; "
    "Yanxiao Han, yanxiao.han@aliyun.com; Yanshuo Han, yanshuohan@dlut.edu.cn. "
    "Yanshuo Han is the submission-system corresponding author.")
letter.core_properties.title = "CSBJ-D-26-00876 revision cover letter"
letter.core_properties.author = "Yanshuo Han"
letter_path = COVER / "CSBJ_cover_letter.docx"
letter.save(letter_path)
(COVER / "CSBJ_cover_letter.txt").write_text(
    "\n\n".join(p.text for p in letter.paragraphs) + "\n", encoding="utf-8")

# Script 31 builds the response from Markdown.
profile = (TMP / "libreoffice_profile").as_uri()
subprocess.run([str(OFFICE), "-env:UserInstallation=" + profile, "--headless",
                "--convert-to", "pdf", "--outdir", str(letter_path.parent), str(letter_path)],
               check=True, capture_output=True, creationflags=subprocess.CREATE_NO_WINDOW)
print("Rendered", letter_path.with_suffix(".pdf").name, flush=True)

# Highlight additions and rewritten blocks directly in the clean OOXML.
# The journal's marked-version file remains fully editable.
from zipfile import ZipFile
from lxml import etree
from difflib import SequenceMatcher
import re
import unicodedata

MARKED = MAN / "CSBJ_SGCB_manuscript_marked.docx"
SHADE = "E8F1FB"
NS = {
    "w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main",
    "m": "http://schemas.openxmlformats.org/officeDocument/2006/math",
    "a": "http://schemas.openxmlformats.org/drawingml/2006/main",
    "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
}
W = "{" + NS["w"] + "}"
R = "{" + NS["r"] + "}"

def read_document(path):
    with ZipFile(path) as archive:
        entries = [(entry, archive.read(entry.filename)) for entry in archive.infolist()]
    parts = {entry.filename: data for entry, data in entries}
    document = etree.fromstring(parts["word/document.xml"])
    relations = etree.fromstring(parts["word/_rels/document.xml.rels"])
    targets = {item.get("Id"): item.get("Target") for item in relations}
    return entries, parts, document, targets

def block_key(block, parts, targets):
    text = "".join(block.xpath(".//w:t/text()|.//m:t/text()", namespaces=NS))
    text = re.sub(r"\s+", " ", unicodedata.normalize("NFKC", text)).strip()
    pictures = tuple(parts["word/" + targets[node.get(R + "embed")]]
        for node in block.xpath(".//a:blip[@r:embed]", namespaces=NS))
    return etree.QName(block).localname, text, pictures

_, original_parts, original, original_targets = read_document(
    ORIGINAL / "manuscript" / "CSBJ_SGCB_manuscript.docx")
entries, clean_parts, clean, clean_targets = read_document(MAN / "CSBJ_SGCB_manuscript.docx")
original_blocks = list(original.find("w:body", NS))[:-1]
clean_blocks = list(clean.find("w:body", NS))[:-1]
original_keys = [block_key(block, original_parts, original_targets) for block in original_blocks]
clean_keys = [block_key(block, clean_parts, clean_targets) for block in clean_blocks]
unchanged = set()
for kind, _, _, start, end in SequenceMatcher(None, original_keys, clean_keys, autojunk=False).get_opcodes():
    if kind == "equal":
        unchanged.update(range(start, end))
changed = [index for index in range(len(clean_blocks)) if index not in unchanged]
for index in changed:
    block = clean_blocks[index]
    paragraphs = [block] if block.tag == W + "p" else list(block.iter(W + "p"))
    for paragraph in paragraphs:
        properties = paragraph.find("w:pPr", NS)
        if properties is None:
            properties = etree.Element(W + "pPr")
            paragraph.insert(0, properties)
        shading = properties.find("w:shd", NS)
        if shading is None:
            shading = etree.Element(W + "shd")
            # CT_PPr places paragraph shading after borders, before tabs/spacing.
            before_shading = {"pStyle", "keepNext", "keepLines", "pageBreakBefore",
                "framePr", "widowControl", "numPr", "suppressLineNumbers", "pBdr"}
            position = next((i for i, child in enumerate(properties)
                if etree.QName(child).localname not in before_shading), len(properties))
            properties.insert(position, shading)
        shading.set(W + "val", "clear")
        shading.set(W + "color", "auto")
        shading.set(W + "fill", SHADE)
with ZipFile(MARKED, "w") as output:
    for entry, data in entries:
        output.writestr(entry, etree.tostring(clean, xml_declaration=True,
            encoding="UTF-8", standalone=True) if entry.filename == "word/document.xml" else data)
_, marked_parts, marked, marked_targets = read_document(MARKED)
marked_blocks = list(marked.find("w:body", NS))[:-1]
assert clean_keys == [block_key(block, marked_parts, marked_targets) for block in marked_blocks]
assert all(marked_parts[name] == data for name, data in clean_parts.items()
           if name != "word/document.xml")
# Compare every editable formula and all block contents, beyond the displayed text.
source_document = etree.fromstring(clean_parts["word/document.xml"])
assert [etree.tostring(node, method="c14n") for node in source_document.findall(".//m:oMath", NS)] == [
    etree.tostring(node, method="c14n") for node in marked.findall(".//m:oMath", NS)]
print("Marked revised blocks", len(changed), "of", len(clean_blocks), flush=True)
subprocess.run([str(OFFICE), "-env:UserInstallation=" + profile, "--headless",
                "--convert-to", "pdf", "--outdir", str(MARKED.parent), str(MARKED)],
               check=True, capture_output=True, creationflags=subprocess.CREATE_NO_WINDOW)
print("Prepared cover letter and highlighted marked manuscript from the clean source.", flush=True)


# Assemble the current journal files and their separate reproducibility archive.
from shutil import copy2
from zipfile import ZIP_DEFLATED
from itertools import product

SUBMISSION = ROOT / "submission"
publication_sources = [
    MAN / "CSBJ_SGCB_manuscript.docx", MAN / "CSBJ_SGCB_manuscript.pdf",
    MAN / "CSBJ_SGCB_manuscript_marked.docx", MAN / "CSBJ_SGCB_manuscript_marked.pdf",
    MAN / "CSBJ_SGCB_supplementary.docx", MAN / "CSBJ_SGCB_supplementary.pdf",
    MAN / "Supplementary_Tables.xlsx",
    ROOT / "response" / "Response_to_Editor_and_Reviewers.docx",
    ROOT / "response" / "Response_to_Editor_and_Reviewers.pdf",
    COVER / "CSBJ_cover_letter.docx", COVER / "CSBJ_cover_letter.pdf"]
publication_targets = [
    SUBMISSION / p.name if p.parent == MAN else SUBMISSION / p.parent.name / p.name
    for p in publication_sources]
figure_names = ["Figure_" + str(n) for n in range(1,7)] + [
    "Figure_S" + str(n) for n in range(1,9)]
figure_sources = [MAN / "figures" / (name+"."+extension)
    for name,extension in product(figure_names,("pdf","svg","png"))]
figure_targets = [SUBMISSION / "figures" / p.name for p in figure_sources]
for src,dest in zip(publication_sources+figure_sources,publication_targets+figure_targets):
    dest.parent.mkdir(parents=True,exist_ok=True)
    copy2(src,dest)
package_archive = ROOT / "release" / f"SGCB_{VERSION}.tar.gz"
review_files = publication_targets + figure_targets
review_archive = ROOT / "release" / f"SGCB_CSBJ_revision_{VERSION}_submission.zip"
with ZipFile(review_archive,"w",compression=ZIP_DEFLATED,compresslevel=6,allowZip64=True) as output:
    for path in review_files:
        output.write(path,path.relative_to(ROOT).as_posix())

# Current experiment entries and fixed evidence; no development caches are copied.
experiment_names = [
    "45_run_current_method.R","123_evaluate_current_api.R","217_cv_matched_families.R",
    "08_recount_existing_null_results.R","20_kirc_patient_validation.R",
    "21_external_kirc_validation.R","28_kirc_score_sensitivity_and_examples.R",
    "48_kirc_distribution_stability.R","66_kirc_current_enrichment.R",
    "95_population_biology_comparison.R","80_new_biology_enrichment.R",
    "81_new_biology_examples.R","sle_classical_monocytes_biology.R","kirc_variability_score.R",
    "25_existing_nb_simulation.R","26_pseudobulk_paired_sensitivity.R",
    "103_proteomics_matched_input.R","231_independent_null_comparators.R",
    "15_prepare_seqc_and_tcga.R","16_evaluate_seqc_qpcr.py","29_publication_tables.py"]
code_files = [ROOT/"experiments"/name for name in experiment_names] + [
    ROOT/"visualization"/"30_publication_figures.py",
    ROOT/"visualization"/"82_new_publication_figures.py",
    ROOT/"visualization"/"overview.py",
    ROOT/"visualization"/"_fonts"/"Inter-400.ttf",
    ROOT/"visualization"/"_fonts"/"Inter-600.ttf",
    ROOT/"visualization"/"_fonts"/"Montserrat-700.ttf",
    ROOT/"writing"/"31_build_documents.py",ROOT/"writing"/"32_render_document_pages.py",
    ROOT/"writing"/"70_prepare_revision_uploads.py",
    ROOT/"development"/"00_current_setup.R",ROOT/"development"/"57_cv_moment_regression.cpp"]
case_labels = ["sle_classical_monocytes_2021","kirc","masld_liver_2023"]
benchmark_labels = [f"{panel}_n{design}" for panel,design in product(
    ("channel","zero"),("10_10","10_50","50_50"))]
null_labels = [f"readcount_null_{n}"+("_streamed" if n==3 else "") for n in range(1,11)]
current_prefixes = [f"sgcb_{VERSION}_"+label for label in benchmark_labels+null_labels+case_labels]
current_prefixes += [f"sgcb_{VERSION}_readcount_null_"]
mean_labels = ["seqc","pseudobulk_kang18","protein_ecoli_lfq","protein_cptac_ups1",
    "protein_sgsds_ratio2","protein_sgsds_ratio2.5"] + [
    f"nb_sim_n{n}_pDE10_ef{effect}" for n,effect in product((3,5,10),("1.5","2.0"))]
results = ROOT / "results"
evidence = set()
for prefix in current_prefixes + ["cv_matched_families_"+VERSION+"_"] + [
        "sgcb_0.5.7_"+label for label in mean_labels] + ["current_masld_liver_2023_"]:
    evidence.update(results.glob(prefix+"*"))
for name in ("30_publication_figures.py","82_new_publication_figures.py"):
    text = (ROOT/"visualization"/name).read_text(encoding="utf-8")
    evidence.update(results/p for p in re.findall(r'R\s*/\s*"([^"]+\.(?:csv|tsv))"',text))
text = (ROOT/"experiments"/"29_publication_tables.py").read_text(encoding="utf-8")
evidence.update(results/p for p in re.findall(r'R\s*/\s*"([^"]+\.(?:csv|tsv))"',text))
evidence.update(results/p for p in [
    "mechanism_example_observations.csv","mechanism_example_summary.csv",
    "taqman_reference_all_assays.csv",
    "population_kirc.rds","covariate_CV2_kirc_discovery_design.rds",
    "population_kirc_covariate_genes.csv","population_kirc_covariate_summary.csv",
    "kirc_discovery_sample_membership.csv"])
for prefix in ("protein_ecoli_lfq","protein_cptac_ups1","protein_sgsds_ratio2","protein_sgsds_ratio2.5"):
    evidence.update(results.glob(prefix+"_matched_preprocessing.rds"))
    evidence.update(results.glob(prefix+"_preprocessing_diagnosis_ratio_mode.rds"))
    evidence.update(results.glob(prefix+"_common_genes.csv"))
native_results = ROOT/"tmp"/"revision_0.5.5_integrated"/"results"
evidence.update(native_results.glob("readcount_independent_null_comparator_*"))

data_files = [ROOT/"data"/name for name in [
    "channel_simulation.rds","zero_calibration.rds","kirc_independent_discovery.rds",
    "kirc_published_purity.csv","ccRCC_exp_log_quantile_normalized.txt",
    "sato_external_clinical.csv","seqc_input.rds","seqc_archived_comparator_scores.csv",
    "GSE56457_SEQC_qPCR_GEOSub.txt.gz"]]
data_files += list((ROOT/"data"/"gtex_readcount_independent_null").glob("split_*.rds"))
data_files += [ROOT/"data"/"tcga_kirc_cache"/name for name in ("tcga_kirc.rds","tcga_kirc_clinical.rds")]
for relative in ("data/new_biology/sle_classical_monocytes_2021",
        "data/new_biology/masld_liver_2023","data/external_methods/DiffVar",
        "data/external_methods/clrDV","data/biological_story_annotation_cache"):
    data_files += [p for p in (ROOT/relative).rglob("*") if p.is_file()]
original_data = ROOT/"data"/"original_inputs"/"benchmark_v2"
data_files += [original_data/"01_bulk_rnaseq"/"data"/"bulk"/"02_simulation"/
    f"sim_n{n}_pDE10_ef{effect}.rds" for n,effect in product((3,5,10),("1.5","2.0"))]
data_files += [original_data/"01_bulk_rnaseq"/"data"/"bulk"/"05_pseudobulk"/"kang18_pseudobulk.rds"]
data_files += [original_data/"03_proteomics"/"data"/relative for relative in (
    "ecoli_lfq_deqms/proteinGroups.txt","sgsds_pxd002370/txt/proteinGroups.txt",
    "sgsds_pxd002370/ratio2.5/txt/proteinGroups.txt","cptac_ups1/proteinGroups.txt")]
data_files += [ROOT/"data"/"original_inputs"/"CSBJ_submission"/"reproducibility"/
    "scripts"/"analysis"/"03_proteomics"/"05_load_proteomics_data.R"]
source_files = [p for p in (ROOT/"software"/"SGCB").rglob("*")
    if p.is_file() and p.suffix.lower() not in {".o",".dll",".so",".a",".d"}]
manuscript_sources = [MAN/"CSBJ_SGCB_manuscript.md",MAN/"CSBJ_SGCB_supplementary.md",
    MAN/"figure_captions.md",MAN/"tables.md",
    ROOT/"response"/"Response_to_Editor_and_Reviewers.md",
    ORIGINAL/"manuscript"/"CSBJ_SGCB_manuscript.docx",
    ORIGINAL/"cover_letter"/"CSBJ_cover_letter.docx"]
reproduction_files = sorted(set(code_files+data_files+source_files+manuscript_sources+
    list(evidence)+[ROOT/"REPRODUCE.md",package_archive]))
reproduction_archive = ROOT/"release"/f"SGCB_CSBJ_revision_{VERSION}_reproducibility.zip"
with ZipFile(reproduction_archive,"w",compression=ZIP_DEFLATED,compresslevel=6,allowZip64=True) as output:
    for path in reproduction_files:
        output.write(path,path.relative_to(ROOT).as_posix())
print("Prepared current submission and reproducibility archives:",review_archive.name,
    reproduction_archive.name,flush=True)
