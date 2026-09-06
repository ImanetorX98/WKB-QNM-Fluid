from __future__ import annotations

import os
import re
from dataclasses import dataclass
from pathlib import Path

from reportlab.lib import colors
from reportlab.lib.colors import HexColor
from reportlab.lib.enums import TA_CENTER, TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import cm, mm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (
    BaseDocTemplate,
    Flowable,
    Frame,
    HRFlowable,
    KeepTogether,
    PageBreak,
    PageTemplate,
    Paragraph,
    Spacer,
    Table,
    TableStyle,
)


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "output" / "pdf"
OUT.mkdir(parents=True, exist_ok=True)
PDF_PATH = OUT / "rapporto_wkb_qnm_madelung.pdf"


NAVY = HexColor("#14263D")
BLUE = HexColor("#246B8E")
TEAL = HexColor("#2A8C82")
PALE = HexColor("#EAF2F4")
GOLD = HexColor("#C9982D")
INK = HexColor("#1D2730")
MID = HexColor("#56636D")
LIGHT = HexColor("#EEF1F3")
WHITE = colors.white


def register_fonts() -> tuple[str, str, str]:
    bundled_fonts = "/Users/iman.rosignoli/.cache/codex-runtimes/codex-primary-runtime/dependencies/native/libreoffice-headless/libreoffice/LibreOfficeDev.app/Contents/Resources/fonts/truetype"
    candidates = [
        (
            f"{bundled_fonts}/DejaVuSans.ttf",
            f"{bundled_fonts}/DejaVuSans-Bold.ttf",
            f"{bundled_fonts}/DejaVuSansMono.ttf",
        ),
        (
            "/System/Library/Fonts/Supplemental/Arial.ttf",
            "/System/Library/Fonts/Supplemental/Arial Bold.ttf",
            "/System/Library/Fonts/Supplemental/Courier New.ttf",
        ),
        (
            "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
            "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
            "/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf",
        ),
    ]
    for regular, bold, mono in candidates:
        if all(Path(p).exists() for p in (regular, bold, mono)):
            pdfmetrics.registerFont(TTFont("ReportRegular", regular))
            pdfmetrics.registerFont(TTFont("ReportBold", bold))
            pdfmetrics.registerFont(TTFont("ReportMono", mono))
            return "ReportRegular", "ReportBold", "ReportMono"
    return "Helvetica", "Helvetica-Bold", "Courier"


FONT, FONT_BOLD, FONT_MONO = register_fonts()


class NumberedCanvasMixin:
    pass


class ReportDocTemplate(BaseDocTemplate):
    def __init__(self, filename: str):
        super().__init__(
            filename,
            pagesize=A4,
            leftMargin=18 * mm,
            rightMargin=18 * mm,
            topMargin=21 * mm,
            bottomMargin=18 * mm,
            title="Interpretazione idrodinamica della WKB per i QNM",
            author="Rapporto di ricerca assistito da Codex",
            subject="WKB, modi quasi-normali, Madelung, exact WKB",
        )
        frame = Frame(
            self.leftMargin,
            self.bottomMargin,
            self.width,
            self.height,
            id="body",
            leftPadding=0,
            rightPadding=0,
            topPadding=0,
            bottomPadding=0,
        )
        self.addPageTemplates(PageTemplate(id="main", frames=[frame], onPage=self._on_page))

    def _on_page(self, canvas, doc):
        if doc.page == 1:
            return
        canvas.saveState()
        canvas.setStrokeColor(HexColor("#D3DADF"))
        canvas.setLineWidth(0.5)
        canvas.line(doc.leftMargin, 14.5 * mm, A4[0] - doc.rightMargin, 14.5 * mm)
        canvas.setFont(FONT, 7.5)
        canvas.setFillColor(MID)
        canvas.drawString(doc.leftMargin, 9.5 * mm, "WKB · QNM · Madelung — rapporto di ricerca")
        canvas.drawRightString(A4[0] - doc.rightMargin, 9.5 * mm, str(doc.page))
        canvas.restoreState()


class VerdictCard(Flowable):
    def __init__(self, width: float, label: str, title: str, body: str, color=TEAL):
        super().__init__()
        self.width = width
        self.height = 31 * mm
        self.label = label
        self.title = title
        self.body = body
        self.color = color

    def draw(self):
        c = self.canv
        c.setFillColor(PALE)
        c.roundRect(0, 0, self.width, self.height, 3 * mm, stroke=0, fill=1)
        c.setFillColor(self.color)
        c.roundRect(0, 0, 6 * mm, self.height, 3 * mm, stroke=0, fill=1)
        c.rect(3 * mm, 0, 3 * mm, self.height, stroke=0, fill=1)
        c.setFillColor(NAVY)
        c.setFont(FONT_BOLD, 8)
        c.drawString(10 * mm, self.height - 8 * mm, self.label.upper())
        c.setFont(FONT_BOLD, 11)
        c.drawString(10 * mm, self.height - 14 * mm, self.title)
        text = c.beginText(10 * mm, self.height - 20 * mm)
        text.setFont(FONT, 8.3)
        text.setFillColor(INK)
        for line in wrap_plain(self.body, 89):
            text.textLine(line)
        c.drawText(text)


class HierarchyDiagram(Flowable):
    def __init__(self, width: float):
        super().__init__()
        self.width = width
        self.height = 60 * mm

    def draw_box(self, c, x, y, w, h, title, body, fill, stroke):
        c.setFillColor(fill)
        c.setStrokeColor(stroke)
        c.setLineWidth(0.8)
        c.roundRect(x, y, w, h, 2.5 * mm, fill=1, stroke=1)
        c.setFillColor(NAVY)
        c.setFont(FONT_BOLD, 8.3)
        c.drawCentredString(x + w / 2, y + h - 6 * mm, title)
        c.setFont(FONT, 7.2)
        c.setFillColor(INK)
        lines = wrap_plain(body, 34)
        ty = y + h - 11 * mm
        for line in lines[:4]:
            c.drawCentredString(x + w / 2, ty, line)
            ty -= 4 * mm

    def arrow(self, c, x1, y1, x2, y2):
        c.setStrokeColor(BLUE)
        c.setFillColor(BLUE)
        c.setLineWidth(1.2)
        c.line(x1, y1, x2, y2)
        import math

        angle = math.atan2(y2 - y1, x2 - x1)
        size = 2.1 * mm
        pts = []
        for delta in (2.55, -2.55):
            pts.append((x2 + size * math.cos(angle + delta), y2 + size * math.sin(angle + delta)))
        p = c.beginPath()
        p.moveTo(x2, y2)
        p.lineTo(*pts[0])
        p.lineTo(*pts[1])
        p.close()
        c.drawPath(p, fill=1, stroke=0)

    def draw(self):
        c = self.canv
        gap = 6 * mm
        w = (self.width - 2 * gap) / 3
        h = 27 * mm
        y_top = 31 * mm
        self.draw_box(c, 0, y_top, w, h, "ODE CANONICA", "ε²ψ″ + qψ = 0\nscelta di gauge", WHITE, BLUE)
        self.draw_box(c, w + gap, y_top, w, h, "RICCATI / MADELUNG", "p² − iεp′ = q\nun solo Q_M[A]", PALE, TEAL)
        self.draw_box(c, 2 * (w + gap), y_top, w, h, "GERARCHIA LOCALE", "uₙ, Q₂ₙ ricorsivi\nnon uniformi ai turning point", WHITE, BLUE)
        self.arrow(c, w, y_top + h / 2, w + gap - 1 * mm, y_top + h / 2)
        self.arrow(c, 2 * w + gap, y_top + h / 2, 2 * (w + gap) - 1 * mm, y_top + h / 2)
        bottom_w = (self.width - gap) / 2
        self.draw_box(c, 0, 0, bottom_w, 23 * mm, "MATCHING DI BARRIERA", "forma normale → Λₙ\ncoefficienti scalari", HexColor("#FFF8E6"), GOLD)
        self.draw_box(c, bottom_w + gap, 0, bottom_w, 23 * mm, "EXACT WKB GLOBALE", "Stokes, Voros symbols\ne periodi quantistici", HexColor("#F0F1FA"), HexColor("#6A6FB0"))
        self.arrow(c, self.width / 2, y_top, bottom_w / 2, 23 * mm)
        self.arrow(c, self.width / 2, y_top, bottom_w + gap + bottom_w / 2, 23 * mm)


def wrap_plain(text: str, width: int) -> list[str]:
    words = text.replace("\n", " \n ").split()
    lines, current = [], ""
    for word in words:
        if word == "\\n":
            lines.append(current)
            current = ""
            continue
        trial = f"{current} {word}".strip()
        if len(trial) > width and current:
            lines.append(current)
            current = word
        else:
            current = trial
    if current:
        lines.append(current)
    return lines


def styles():
    ss = getSampleStyleSheet()
    return {
        "body": ParagraphStyle(
            "Body",
            parent=ss["BodyText"],
            fontName=FONT,
            fontSize=9.25,
            leading=13.1,
            textColor=INK,
            spaceAfter=4.5,
            allowWidows=0,
            allowOrphans=0,
        ),
        "small": ParagraphStyle(
            "Small",
            parent=ss["BodyText"],
            fontName=FONT,
            fontSize=7.4,
            leading=10.2,
            textColor=INK,
        ),
        "caption": ParagraphStyle(
            "Caption",
            parent=ss["BodyText"],
            fontName=FONT,
            fontSize=7.5,
            leading=10,
            textColor=MID,
            alignment=TA_CENTER,
            spaceBefore=2,
            spaceAfter=7,
        ),
        "h1": ParagraphStyle(
            "H1",
            parent=ss["Heading1"],
            fontName=FONT_BOLD,
            fontSize=17,
            leading=20,
            textColor=NAVY,
            spaceBefore=12,
            spaceAfter=7,
            keepWithNext=1,
        ),
        "h2": ParagraphStyle(
            "H2",
            parent=ss["Heading2"],
            fontName=FONT_BOLD,
            fontSize=12.5,
            leading=15,
            textColor=BLUE,
            spaceBefore=9,
            spaceAfter=5,
            keepWithNext=1,
        ),
        "h3": ParagraphStyle(
            "H3",
            parent=ss["Heading3"],
            fontName=FONT_BOLD,
            fontSize=10.3,
            leading=13,
            textColor=TEAL,
            spaceBefore=6,
            spaceAfter=3,
            keepWithNext=1,
        ),
        "eq": ParagraphStyle(
            "Equation",
            parent=ss["BodyText"],
            fontName=FONT_MONO,
            fontSize=8.1,
            leading=11.5,
            textColor=NAVY,
            backColor=HexColor("#F4F7F8"),
            borderColor=HexColor("#DCE5E8"),
            borderWidth=0.6,
            borderPadding=(5, 7, 5, 7),
            spaceBefore=5,
            spaceAfter=7,
            alignment=TA_CENTER,
        ),
        "callout": ParagraphStyle(
            "Callout",
            parent=ss["BodyText"],
            fontName=FONT,
            fontSize=9.2,
            leading=13,
            textColor=NAVY,
            backColor=PALE,
            borderColor=TEAL,
            borderWidth=1,
            borderPadding=(7, 9, 7, 9),
            spaceBefore=6,
            spaceAfter=8,
        ),
        "bullet": ParagraphStyle(
            "Bullet",
            parent=ss["BodyText"],
            fontName=FONT,
            fontSize=9,
            leading=12.6,
            leftIndent=13,
            firstLineIndent=-7,
            bulletIndent=4,
            textColor=INK,
            spaceAfter=2.5,
        ),
        "ref": ParagraphStyle(
            "Reference",
            parent=ss["BodyText"],
            fontName=FONT,
            fontSize=7.4,
            leading=10.2,
            leftIndent=14,
            firstLineIndent=-14,
            textColor=INK,
            spaceAfter=3.2,
        ),
        "toc": ParagraphStyle(
            "TOC",
            parent=ss["BodyText"],
            fontName=FONT,
            fontSize=9.4,
            leading=15,
            textColor=INK,
        ),
    }


S = styles()


def P(text: str, style="body") -> Paragraph:
    return Paragraph(text, S[style])


def Eq(text: str) -> Paragraph:
    return P(text, "eq")


def bullets(items: list[str]):
    return [P("• " + item, "bullet") for item in items]


def section(title: str):
    return [Spacer(1, 2 * mm), HRFlowable(width="100%", thickness=0.7, color=HexColor("#CAD5DA"), spaceAfter=4), P(title, "h1")]


def subsection(title: str):
    return [P(title, "h2")]


def table(data, widths, header=True, font_size=7.4, leading=9.2):
    cooked = []
    for r, row in enumerate(data):
        cooked.append([
            cell if isinstance(cell, Flowable) else Paragraph(str(cell), ParagraphStyle(
                f"cell-{r}", fontName=FONT_BOLD if header and r == 0 else FONT,
                fontSize=font_size, leading=leading,
                textColor=WHITE if header and r == 0 else INK,
            ))
            for cell in row
        ])
    t = Table(cooked, colWidths=widths, repeatRows=1 if header else 0, hAlign="LEFT")
    style = [
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LEFTPADDING", (0, 0), (-1, -1), 5),
        ("RIGHTPADDING", (0, 0), (-1, -1), 5),
        ("TOPPADDING", (0, 0), (-1, -1), 5),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
        ("GRID", (0, 0), (-1, -1), 0.35, HexColor("#CCD6DB")),
    ]
    if header:
        style += [("BACKGROUND", (0, 0), (-1, 0), NAVY)]
        if len(data) > 1:
            style += [("ROWBACKGROUNDS", (0, 1), (-1, -1), [WHITE, HexColor("#F7F9FA")])]
    t.setStyle(TableStyle(style))
    return t


def link(url: str, label: str) -> str:
    return f'<a href="{url}" color="#246B8E"><u>{label}</u></a>'


def build_story():
    story = []

    # Cover
    story += [Spacer(1, 19 * mm)]
    story.append(P("RAPPORTO DI RICERCA", "h3"))
    story.append(Spacer(1, 7 * mm))
    cover_title = ParagraphStyle("CoverTitle", fontName=FONT_BOLD, fontSize=27, leading=31, textColor=NAVY, spaceAfter=10)
    cover_sub = ParagraphStyle("CoverSub", fontName=FONT, fontSize=13, leading=18, textColor=BLUE, spaceAfter=8)
    story.append(Paragraph("Interpretazione idrodinamica<br/>dell’espansione WKB", cover_title))
    story.append(Paragraph("Modi quasi-normali dei buchi neri e formalismo di Madelung", cover_sub))
    story.append(Spacer(1, 8 * mm))
    story.append(HRFlowable(width="28%", thickness=2.2, color=TEAL, hAlign="LEFT"))
    story.append(Spacer(1, 8 * mm))
    story.append(P("<b>Obiettivo.</b> Verificare — senza assumerla — l’esistenza di una gerarchia idrodinamica all-order, reinterpretando le formule WKB già stabilite e distinguendo identità formali, quantizzazione di barriera e contenuto fisico.", "callout"))
    story.append(Spacer(1, 13 * mm))
    meta = [
        ["Stato", "Sintesi critica con formulazione matematica e programma falsificabile"],
        ["Letteratura", "Schutz–Will · Iyer–Will · Konoplya · Matyjasek–Opala · exact WKB"],
        ["Ricerca chiusa", "23 agosto 2026"],
    ]
    story.append(table(meta, [31 * mm, 126 * mm], header=False, font_size=8.2, leading=11))
    story.append(Spacer(1, 17 * mm))
    story.append(P("Conclusione in una riga", "h3"))
    story.append(Paragraph("La gerarchia esiste come espansione di <b>un unico funzionale di Madelung</b>; non emerge invece una famiglia canonica di potenziali locali identificabile ordine per ordine con le correzioni QNM di barriera.", ParagraphStyle("CoverVerdict", fontName=FONT_BOLD, fontSize=12, leading=17, textColor=NAVY)))
    story.append(PageBreak())

    # TOC
    story += section("Mappa del rapporto")
    toc_rows = [
        ("01", "Domande e livelli dell’ipotesi"),
        ("02", "Convenzioni e dominio di validità"),
        ("03", "WKB/QNM: struttura usata"),
        ("04", "Chiusura locale di Madelung a tutti gli ordini"),
        ("05", "Frequenze complesse e legge di continuità"),
        ("06", "Ostruzioni alla mappa locale ↔ Λₙ"),
        ("07", "Significato delle correzioni dispersive"),
        ("08", "Revisione della letteratura e lacune"),
        ("09", "Programma falsificabile"),
        ("10", "Applicazioni e conclusione"),
        ("11", "Registro delle affermazioni e bibliografia"),
    ]
    toc_data = [[P(f"<b>{n}</b>", "toc"), P(title, "toc")] for n, title in toc_rows]
    toc_table = Table(toc_data, colWidths=[16 * mm, 142 * mm], hAlign="LEFT")
    toc_table.setStyle(TableStyle([
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("LINEBELOW", (0, 0), (-1, -1), 0.35, HexColor("#D7DEE2")),
        ("TOPPADDING", (0, 0), (-1, -1), 6),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 6),
        ("TEXTCOLOR", (0, 0), (0, -1), TEAL),
    ]))
    story.append(toc_table)
    story.append(Spacer(1, 9 * mm))
    story.append(P("<b>Legenda epistemica.</b> “Confermato” indica un’identità o un risultato sostenuto direttamente dalle fonti. “Dedotto” indica un argomento matematico esplicitato nel rapporto. “Proposto” indica una direzione di ricerca ancora da validare.", "callout"))
    story.append(PageBreak())

    # Executive summary
    story += section("Sintesi esecutiva")
    story.append(VerdictCard(159 * mm, "Confermato", "Chiusura formale locale", "Un solo funzionale Q_M[A] genera, per espansione auto-consistente, tutti i termini pari della WKB.", TEAL))
    story.append(Spacer(1, 4 * mm))
    story.append(VerdictCard(159 * mm, "Non confermato", "Potenziali indipendenti per ordine", "I coefficienti locali non sono canonici né identificabili uno-a-uno con i Λₙ del matching di barriera.", GOLD))
    story.append(Spacer(1, 4 * mm))
    story.append(VerdictCard(159 * mm, "Aperto", "Interpretazione di fluido fisico", "Per frequenze complesse serve un fluido con sorgente oppure una complessificazione senza densità positiva ordinaria.", HexColor("#6A6FB0")))
    story.append(Spacer(1, 7 * mm))
    story.append(P("Il risultato centrale è più preciso della congettura iniziale: la gerarchia locale esiste, ma come serie derivata da <b>un’unica pressione quantistica</b>. Il candidato più robusto per una formulazione globale e invariante non è una collezione di potenziali locali: sono il momento quantistico come 1-forma, i periodi e i dati di Stokes dell’exact WKB."))
    story.append(Spacer(1, 5 * mm))
    story.append(HierarchyDiagram(159 * mm))
    story.append(P("Figura 1 — Relazioni e biforcazione concettuale fra gerarchia locale, matching di barriera ed exact WKB globale.", "caption"))

    # 1
    story += section("1. Domande di ricerca e livelli dell’ipotesi")
    story.append(P("La domanda iniziale mescola tre enunciati diversi. Separarli evita di trasformare un’identità algebrica utile in una dichiarazione fisica non dimostrata."))
    hdata = [
        ["Ipotesi", "Contenuto", "Esito"],
        ["H1 · chiusura formale", "La serie WKB completa si esprime con un funzionale di tipo Madelung e una ricorsione compatta.", "Confermata localmente, lontano da zeri e turning point."],
        ["H2 · potenziali per ordine", "Ogni correzione WKB/QNM è un potenziale locale autonomo e canonico.", "Non confermata; la ricostruzione è in generale non unica."],
        ["H3 · fluido fisico", "La gerarchia descrive correzioni dispersive di un fluido reale associato ai QNM.", "Aperta ma vincolata: serve una sorgente o una complessificazione."],
    ]
    story.append(table(hdata, [32 * mm, 72 * mm, 55 * mm], font_size=7.7, leading=10.2))

    # 2
    story += section("2. Convenzioni e dominio di validità")
    story.append(P("Consideriamo la forma canonica senza derivata prima"))
    story.append(Eq("ε² ψ″(x) + q(x,ω) ψ(x) = 0,          q = ω² − V(x)  nel caso base."))
    story.append(P("Il tempo è proporzionale a exp(−iωt); per un QNM stabile Im ω &lt; 0. Il parametro ε è formale e viene posto uguale a uno alla fine. La letteratura exact-WKB recente nota che il modo di introdurlo non è unico: un’interpretazione fisica dei coefficienti deve fissare coordinata tortoise, variabile master, normalizzazione e scala semiclassica."))
    story += bullets([
        "Dominio locale semplicemente connesso con q ≠ 0 e un ramo fissato di √q.",
        "Assenza di zeri di ψ nel dominio della trasformazione di Riccati.",
        "Turning point, linee di Stokes e condizioni al contorno globali richiedono uniformizzazione o continuazione analitica.",
    ])

    # 3
    story += section("3. WKB di barriera per QNM: la struttura reinterpretata")
    story.append(P("Schutz–Will introducono l’approccio semianalitico WKB; Iyer–Will ne costruiscono la sistematica ad ordine superiore mediante matching uniforme al massimo; Konoplya estende il calcolo al sesto ordine; Matyjasek–Opala al tredicesimo con Padé. Sviluppi successivi arrivano a ordini molto elevati e a formulazioni exact-WKB."))
    story.append(Eq("i Q₀ / √(2 Q₀″) − Σₖ₌₂ᴺ Λₖ = n + 1/2"))
    story.append(P("Qui Q = ω² − V; il pedice 0 indica il massimo della barriera. Ogni Λₖ è una combinazione del getto Q₀⁽³⁾,…,Q₀⁽²ᵏ⁾ e di α = n + 1/2."))
    story += bullets([
        "I Λₖ sono coefficienti scalari di matching/quantizzazione, non campi locali.",
        "La loro dipendenza dall’overtone non è quella di un potenziale di fondo universale.",
        "La serie è asintotica; ordine ottimo, Padé o Borel–Padé sono parte del problema, non dettagli numerici accessori.",
    ])

    # 4
    story += section("4. Proposizione: chiusura locale di Madelung a tutti gli ordini")
    story += subsection("4.1 Equazione di Riccati e ricorsione WKB")
    story.append(P("Definiamo il momento quantistico p = −iε ψ′/ψ. Dove ψ non si annulla, l’ODE è esattamente equivalente a"))
    story.append(Eq("p² − i ε p′ = q.                                                     (4.1)"))
    story.append(P("Con p = Σₙ≥₀ εⁿpₙ si ritrova la ricorsione standard:"))
    story.append(Eq("p₀² = q,     pₙ = [ i pₙ₋₁′ − Σₖ₌₁ⁿ⁻¹ pₖ pₙ₋ₖ ] / (2p₀).          (4.2)"))
    story += subsection("4.2 Separazione pari/dispari")
    story.append(P("Scriviamo p = u + v, con u pari e v dispari in ε. La parte dispari fissa v; la parte pari chiude l’equazione per u:"))
    story.append(Eq("v = (iε/2) u′/u,       ψ = C u^(-1/2) exp[(i/ε)∫u dx].             (4.3)"))
    story.append(Eq("q = u² + (ε²/2)u″/u − (3ε²/4)(u′/u)².                             (4.4)"))
    story.append(P("Ponendo A = u^(-1/2) si ottiene il punto chiave:"))
    story.append(Eq("q = u² + Q_M[A],              Q_M[A] = −ε² A″/A.                  (4.5)"))
    story.append(P("<b>Risultato.</b> Tutta la gerarchia WKB locale è generata da un solo funzionale di Madelung valutato sull’ampiezza completa e dipendente da ε. Non serve postulare un nuovo funzionale a ogni ordine.", "callout"))
    story += subsection("4.3 Ricorsione compatta della gerarchia efficace")
    story.append(Eq("u = Σⱼ≥₀ ε²ʲuⱼ,      u₀ = √q,      F[u] = ¾(u′/u)² − ½u″/u."))
    story.append(P("Dalla relazione u² = q + ε²F[u], per n ≥ 1:"))
    story.append(Eq("uₙ = { [ε²ⁿ⁻²]F[u] − Σⱼ₌₁ⁿ⁻¹uⱼuₙ₋ⱼ } / (2u₀).                   (4.6)"))
    story.append(Eq("Q_eff = Σₙ≥₁ ε²ⁿQ₂ₙ,       Q₂ₙ = −[ε²ⁿ⁻²]F[u].                    (4.7)"))
    story.append(P("Il primo livello, utile anche come controllo simbolico, è"))
    story.append(Eq("Q₂ = q″/(4q) − 5(q′)²/(16q²),\n u₁ = 5(q′)²/(32 q^(5/2)) − q″/(8 q^(3/2)).                           (4.8)"))
    story.append(P("La (4.6) è la formulazione generale richiesta: compatta, algoritmica e verificabile ordine per ordine. Tuttavia i coefficienti diventano singolari quando q → 0; ciò non è un piccolo difetto, ma il segnale che la rappresentazione locale non è uniforme."))
    story += subsection("4.4 Portata logica")
    story += bullets([
        "Dimostra una chiusura Madelung all-order della WKB locale.",
        "Non rende Q₂ₙ osservabili indipendenti o invarianti di coordinata.",
        "Non prova l’identità Q₂ₙ ↔ Λₙ.",
        "Non introduce automaticamente nuovi termini microscopici nella dispersione.",
    ])

    # 5
    story += section("5. Frequenze complesse: fluido con sorgente o fluido complesso")
    story.append(P("Imponiamo ora ampiezza e fase reali, ψ = √ρ exp(iθ/ε), e scriviamo q = q_R + iq_I. La separazione reale/immaginaria dà esattamente"))
    story.append(Eq("(θ′)² − ε²(√ρ)″/√ρ = q_R,                                          (5.1)\n(ρθ′)′ = −(q_I/ε)ρ.                                                    (5.2)"))
    story.append(P("La prima è Hamilton–Jacobi con potenziale quantistico; la seconda è una continuità con sorgente o pozzo. Per q = ω² − V e V reale, q_I = Im(ω²)."))
    story.append(P("Con Im ω &lt; 0, le condizioni QNM outgoing/ingoing producono inoltre autofunzioni divergenti alle estremità del dominio tortoise. Non sono stati L². Una densità positiva globale va quindi intesa come oggetto locale o come variabile di un sistema aperto."))
    data = [
        ["Scelta", "Vantaggio", "Costo concettuale"],
        ["ρ, θ reali", "Densità e fase leggibili idrodinamicamente", "Continuità con sorgente; non conservativa"],
        ["u, A complessi", "Compatibilità diretta con exact WKB e condizioni QNM", "Nessuna densità positiva/velocità reale ordinaria"],
    ]
    story.append(table(data, [37 * mm, 58 * mm, 64 * mm], font_size=7.7, leading=10.2))

    # 6
    story += section("6. Ostruzioni a una mappa canonica Q₂ₙ(x) ↔ Λₙ")
    story += subsection("6.1 Non uniformità al massimo della barriera")
    story.append(P("I Q₂ₙ locali contengono inverse potenze di q. Per i QNM a basso overtone i turning point coalescono vicino al massimo e q(x₀) è piccolo nella scalatura semiclassica. Proprio nella regione decisiva per lo spettro, la gerarchia locale diverge."))
    story.append(P("Iyer–Will espandono il potenziale al massimo, riscalano x − x₀ ∼ √ε e riducono il problema a una forma normale di cilindro parabolico perturbata. I Λₙ appartengono a questa procedura uniforme, non sono semplicemente valori dei Q₂ₙ locali."))
    story += subsection("6.2 Non unicità inversa")
    story.append(P("A ordine N, la quantizzazione usa un getto finito del potenziale in x₀. Infinite funzioni δV_N(x) possono condividere quel getto e dare la stessa correzione scalare. Dato Λ_N, il potenziale efficace locale è quindi sottodeterminato senza un ansatz o un gauge aggiuntivo."))
    story += subsection("6.3 Trasformazioni di Liouville")
    story.append(P("Per x = x(z) e φ(z) = ψ(x(z))/√x′(z), la forma canonica è preservata con"))
    story.append(Eq("q̃(z) = (x′)² q(x(z)) + (ε²/2){x,z}.                                (6.1)"))
    story.append(P("La derivata schwarziana può spostare un termine di ordine ε² fra potenziale di fondo e correzione quantistica. Un Q₂ₙ locale non è quindi, da solo, un oggetto geometrico invariante."))
    story.append(P("<b>Deduzione di non unicità.</b> La combinazione di non uniformità, getto finito e termine schwarziano impedisce una mappa canonica Λₙ → Q₂ₙ(x) senza ulteriori convenzioni. Questa è una deduzione del rapporto, non un teorema attribuito agli autori WKB.", "callout"))
    story += subsection("6.4 Oggetti globali più robusti")
    story.append(P("L’exact WKB organizza la soluzione tramite il momento quantistico, le linee di Stokes, i Voros symbols e i periodi su cicli nel piano complesso. Questi oggetti incorporano la continuazione globale e sono candidati più naturali per una formulazione relativamente invariante."))

    # 7
    story += section("7. Correzioni dispersive: una distinzione necessaria")
    story.append(P("Nella Schrödinger standard, il termine di pressione quantistica è già un gradiente di ordine due. La sua iterazione nella soluzione WKB produce derivate sempre più alte di q: questa è una <b>gerarchia asintotica della soluzione</b>."))
    story.append(P("Una vera gerarchia costitutiva di gradienti superiori — ad esempio in modelli non locali di fase o in EFT con operatori ad alte derivate — modifica invece l’equazione dinamica e la relazione di dispersione. I due livelli possono essere confrontati, non identificati senza derivazione."))
    story.append(P("Formulazione ammessa: i Q₂ₙ organizzano correzioni di gradiente/dispersive <i>efficaci della soluzione</i>. Formulazione non ancora giustificata: i Q₂ₙ sono nuovi coefficienti microscopici di dispersione del mezzo.", "callout"))

    # 8
    story += section("8. Revisione della letteratura e matrice delle lacune")
    story += subsection("8.1 Quattro filoni che si sfiorano")
    fields = [
        ["Filone", "Risultato consolidato", "Lacuna rispetto al progetto"],
        ["WKB/QNM", "Ordini elevati, matching uniforme, Padé/Borel, codici", "Nessuna lettura Madelung sistematica dei coefficienti"],
        ["Madelung / QHJ", "Riccati, potenziale quantistico, pari/dispari exact-WKB", "Applicazione esplicita ai QNM gravitazionali di barriera"],
        ["Analog gravity", "QNM di BEC e dispersione Bogoliubov fisica", "Non è una reinterpretazione dei Λₙ gravitazionali"],
        ["Gradienti superiori", "Gerarchie costitutive in modelli non locali", "Equivalenza con la gerarchia WKB non dimostrata"],
    ]
    story.append(table(fields, [31 * mm, 64 * mm, 64 * mm], font_size=7.3, leading=9.7))
    story += subsection("8.2 Ricerca d’intersezione")
    story.append(P("Sono state cercate combinazioni fra “Madelung”, “Bohm quantum potential”, “quasinormal mode”, “black hole”, “higher-order WKB”, “Schutz–Will” e “Iyer–Will” in motori bibliografici, archivi preprint e catene di citazioni. Sono stati controllati i lavori exact-WKB e analog-gravity più vicini."))
    story.append(P("<b>Nessun precedente diretto è stato identificato al 23 agosto 2026.</b> La formulazione è intenzionalmente qualificata: indicizzazione incompleta, terminologie diverse e letteratura non digitalizzata impediscono di trasformare un esito bibliografico in prova di inesistenza.", "callout"))
    gap = [
        ["Affermazione", "Status", "Confidenza"],
        ["WKB QNM ad alto ordine disponibile", "Consolidato", "Alta"],
        ["Chiusura Riccati/Madelung all-order", "Identità locale", "Alta"],
        ["Mappa canonica Λₙ ↔ potenziale locale", "Non dimostrata / non unica", "Alta"],
        ["Precedente diretto Madelung–WKB–QNM", "Non identificato", "Media"],
        ["Nuova idrodinamica fisica dei QNM", "Ipotesi aperta", "Bassa oggi"],
    ]
    story.append(table(gap, [78 * mm, 50 * mm, 31 * mm], font_size=7.6, leading=10))

    # 9
    story += section("9. Programma di ricerca falsificabile")
    phases = [
        ["Fase", "Azioni", "Output / test"],
        ["A · Gauge", "Fissare metrica, tortoise, variabile master, ε e trasformazioni ammesse", "Definizione riproducibile"],
        ["B · Gerarchia locale", "Generare uₙ e Q₂ₙ; verificare Riccati; studiare crescita e singolarità", "Algoritmo simbolico + diagnostica"],
        ["C · Uniform.", "Forma normale parabolico-cilindrica; confronto con Λ₂, Λ₃, …", "Ponte oppure confutazione"],
        ["D · Complessità", "Confrontare fluido con sorgente e chiusura complessa; test Darboux", "Criterio di significato fisico"],
        ["E · Benchmark", "Pöschl–Teller, Schwarzschild, RN, doppia barriera; Padé/Borel", "Errori contro spettro di riferimento"],
    ]
    story.append(table(phases, [25 * mm, 88 * mm, 46 * mm], font_size=7.25, leading=9.7))
    story += subsection("Criteri di confutazione della versione forte")
    story += bullets([
        "I potenziali cambiano sotto trasformazioni ammesse mentre lo spettro resta invariato.",
        "Dopo gauge fissato e uniformizzazione non si riproducono almeno Λ₂ e Λ₃.",
        "Non si definisce densità/flusso coerente per frequenze complesse.",
        "Nel benchmark esatto la gerarchia non organizza né migliora le approssimazioni rispetto alla WKB ordinaria.",
    ])
    story += subsection("Benchmark prioritari")
    story += bullets([
        "Pöschl–Teller: spettro QNM esatto e barriera singola controllata.",
        "Regge–Wheeler/Zerilli: test di isospectralità e dipendenza dalla variabile master.",
        "Reissner–Nordström: confronto con risultati Borel-sommati esistenti.",
        "Doppia barriera o geometria dirty: separazione fra informazione locale e globale.",
    ])

    # 10
    story += section("10. Applicazioni plausibili e conclusione")
    apps = [
        ["Priorità", "Applicazione", "Valutazione"],
        ["1", "Organizzazione simbolica e diagnostica locale |Q_M|/|q|", "Alta plausibilità; non è un bound d’errore"],
        ["2", "Periodi quantistici e risommazione invariante", "Alta plausibilità matematica"],
        ["3", "Sensibilità dello spettro alla geometria/EFT della barriera", "Plausibilità media"],
        ["4", "Test di isospectralità e scelta della variabile master", "Plausibilità media"],
        ["5", "Ponte con analog gravity e BEC", "Esplorativo; non prova il caso gravitazionale"],
        ["6", "Nuova idrodinamica fisica dei QNM", "Bassa allo stato attuale; serve una teoria del mezzo"],
    ]
    story.append(table(apps, [17 * mm, 77 * mm, 65 * mm], font_size=7.4, leading=9.9))
    story.append(Spacer(1, 5 * mm))
    story.append(P("La congettura contiene un nucleo matematico solido: la WKB locale completa ammette una chiusura Madelung esatta e la famiglia Q₂ₙ è generata ricorsivamente da un unico potenziale quantistico. Questo è il risultato positivo."))
    story.append(P("Il risultato negativo è altrettanto importante: non è giustificata una famiglia fisica e canonica di potenziali identificabile con i Λₙ. La regione spettralmente decisiva richiede uniformizzazione; i coefficienti locali dipendono dal gauge; le frequenze QNM rendono il flusso aperto o complesso."))
    story.append(P("<b>Direzione raccomandata.</b> Usare la gerarchia locale come strumento organizzativo; cercare il contenuto invariante nei periodi e nei dati di Stokes dell’exact WKB, oppure in un sistema idrodinamico aperto dichiarato esplicitamente. La proposta resta utile sia se il ponte con Λₙ funziona, sia se i test lo confutano.", "callout"))

    # 11 ledger
    story += section("11. Registro delle affermazioni")
    ledger = [
        ["Affermazione", "Tipo", "Base", "Conf."],
        ["Chiusura WKB nel singolo Q_M", "Identità", "Riccati + pari/dispari", "Alta"],
        ["Nessuna mappa canonica Λₙ → Q₂ₙ(x)", "Deduzione", "Turning point + getto + Schwarziana", "Alta"],
        ["Continuità con sorgente per q complesso", "Identità", "Separazione reale/immaginaria", "Alta"],
        ["Q₂ₙ non implica nuova micro-dispersione", "Distinzione", "Ordine ODE vs ordine asintotico", "Alta"],
        ["Nessun precedente diretto identificato", "Esito ricerca", "Query e catene di citazioni", "Media"],
        ["Periodi quantistici come candidati invarianti", "Proposta", "Exact WKB + Liouville", "Medio-alta"],
        ["Nuova idrodinamica fisica QNM", "Ipotesi", "Non dimostrata", "Bassa"],
    ]
    story.append(table(ledger, [58 * mm, 27 * mm, 51 * mm, 23 * mm], font_size=6.9, leading=9.2))

    # References
    story += section("Bibliografia essenziale")
    refs = [
        (1, "Schutz, B. F.; Will, C. M. (1985), “Black Hole Normal Modes: A Semianalytic Approach”, ApJL 291 L33.", "https://doi.org/10.1086/184453"),
        (2, "Iyer, S.; Will, C. M. (1987), “Black-hole normal modes: A WKB approach. I”, PRD 35, 3621.", "https://doi.org/10.1103/PhysRevD.35.3621"),
        (3, "Iyer, S. (1987), “Black-hole normal modes: A WKB approach. II”, PRD 35, 3632.", "https://doi.org/10.1103/PhysRevD.35.3632"),
        (4, "Konoplya, R. A. (2003), higher-order WKB for D-dimensional Schwarzschild QNM, PRD 68, 024018.", "https://arxiv.org/abs/gr-qc/0303052"),
        (5, "Matyjasek, J.; Opala, M. (2017), improved semianalytic approach, 13th-order WKB and Padé, PRD 96, 024011.", "https://arxiv.org/abs/1704.00361"),
        (6, "Konoplya, R. A.; Zhidenko, A.; Zinhailo, A. F. (2019), higher-order WKB recipes, CQG 36.", "https://arxiv.org/abs/1904.10333"),
        (7, "Matyjasek, J.; Telecka, M. (2019), Padé summation of very high-order WKB terms, PRD 100, 124006.", "https://arxiv.org/abs/1908.09389"),
        (8, "Hatsuda, Y. (2020), QNM and Borel summation, PRD 101, 024008.", "https://arxiv.org/abs/1906.07232"),
        (9, "Hatsuda, Y.; Kimura, M. (2021), “Spectral Problems for Quasinormal Modes of Black Holes”.", "https://arxiv.org/abs/2111.15197"),
        (10, "Miyachi, T.; Namba, R.; Omiya, H.; Oshita, N. (2025), exact-WKB path for black-hole QNM.", "https://arxiv.org/abs/2503.17245"),
        (11, "Hatsuda, Y.; Shiga, Y. (2026), exact WKB and quantum periods for extremal black-hole QNM.", "https://arxiv.org/abs/2605.01321"),
        (12, "Konoplya, R. A.; Matyjasek, J.; Zhidenko, A. (2026), efficient higher-order WKB code.", "https://arxiv.org/abs/2603.12466"),
        (13, "Madelung, E. (1927), “Quantentheorie in hydrodynamischer Form”, Z. Phys. 40, 322–326.", "https://doi.org/10.1007/BF01400372"),
        (14, "Takabayasi, T. (1952), quantum mechanics associated with classical pictures, PTP 8, 143.", "https://doi.org/10.1143/ptp/8.2.143"),
        (15, "Carles, R.; Danchin, R.; Saut, J.-C. (2012), “Madelung, Gross–Pitaevskii and Korteweg”, Nonlinearity 25.", "https://arxiv.org/abs/1111.4670"),
        (16, "Türe, M.; Ünsal, M. (2025), Quantum Hamilton–Jacobi theory and exact WKB, PRD 111, 105010.", "https://arxiv.org/abs/2406.07829"),
        (17, "Silva, H. O. et al. (2024), phase-amplitude/Riccati methods for black-hole ringdown, PRD 110, 024042.", "https://arxiv.org/abs/2404.11110"),
        (18, "Barceló, C. et al. (2007), QNM in BEC acoustic black holes, PRD 75, 084024.", "https://arxiv.org/abs/gr-qc/0701173"),
        (19, "Daghigh, R. G.; Green, M. D. (2015), high-overtone QNM of analog black holes.", "https://arxiv.org/abs/1411.7066"),
        (20, "Mauri, R. (2021), non-local phase-field model of Bohm’s quantum potential, Found. Phys. 51.", "https://doi.org/10.1007/s10701-021-00454-9"),
        (21, "NIST DLMF §1.13, Liouville transformation.", "https://dlmf.nist.gov/1.13"),
    ]
    for n, citation, url in refs:
        story.append(P(f"[{n}] {citation} {link(url, 'link')}", "ref"))

    story += section("Nota metodologica")
    story.append(P("La revisione privilegia articoli originali, DOI, arXiv e documentazione matematica primaria. Le affermazioni negative sono formulate come “non identificato”, non come “non esiste”. Le deduzioni originali sono marcate. Non è assunta validità uniforme della WKB locale presso i turning point."))
    story.append(P("Le formule sono espresse in convenzioni adimensionali. Fattori di massa e ℏ possono essere ripristinati in base all’equazione master scelta; farlo prima di fissare il gauge oscurerebbe, anziché chiarire, il contenuto invariante."))
    story.append(Spacer(1, 8 * mm))
    story.append(HRFlowable(width="100%", thickness=1.2, color=TEAL))
    story.append(P("Fine del rapporto · ricerca chiusa il 23 agosto 2026", "caption"))
    return story


def main():
    doc = ReportDocTemplate(str(PDF_PATH))
    doc.build(build_story())
    print(PDF_PATH)


if __name__ == "__main__":
    main()
