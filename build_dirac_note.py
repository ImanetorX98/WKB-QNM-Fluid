from __future__ import annotations

from pathlib import Path
from reportlab.lib import colors
from reportlab.lib.colors import HexColor
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import mm
from reportlab.platypus import (
    BaseDocTemplate,
    Frame,
    HRFlowable,
    Image,
    KeepTogether,
    PageBreak,
    PageTemplate,
    Paragraph,
    Spacer,
)

from build_report import (
    BLUE,
    FONT,
    FONT_BOLD,
    GOLD,
    INK,
    MID,
    NAVY,
    PALE,
    TEAL,
    WHITE,
    Eq,
    P,
    VerdictCard,
    bullets,
    link,
    section,
    subsection,
    table,
)


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "output" / "pdf"
TMP = ROOT / "tmp" / "pdfs"
OUT.mkdir(parents=True, exist_ok=True)
TMP.mkdir(parents=True, exist_ok=True)

PDF_PATH = OUT / "nota_dirac_wkb_idrodinamica.pdf"
PLOT_PATH = TMP / "dirac_partner_potentials.png"


class DiracNoteDocTemplate(BaseDocTemplate):
    def __init__(self, filename: str):
        super().__init__(
            filename,
            pagesize=A4,
            leftMargin=18 * mm,
            rightMargin=18 * mm,
            topMargin=21 * mm,
            bottomMargin=18 * mm,
            title="Dirac WKB idrodinamica su Schwarzschild e Kerr",
            author="Nota di calcolo assistita da Codex",
            subject="Dirac, WKB, Madelung, Schwarzschild, Kerr, QNM",
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
        self.addPageTemplates(
            PageTemplate(id="main", frames=[frame], onPage=self._on_page)
        )

    def _on_page(self, canvas, doc):
        if doc.page == 1:
            return
        canvas.saveState()
        canvas.setStrokeColor(HexColor("#D3DADF"))
        canvas.setLineWidth(0.5)
        canvas.line(doc.leftMargin, 14.5 * mm, A4[0] - doc.rightMargin, 14.5 * mm)
        canvas.setFont(FONT, 7.4)
        canvas.setFillColor(MID)
        canvas.drawString(
            doc.leftMargin, 9.5 * mm, "Dirac WKB idrodinamica - nota di calcolo"
        )
        canvas.drawRightString(
            A4[0] - doc.rightMargin, 9.5 * mm, str(doc.page)
        )
        canvas.restoreState()


def make_potential_plot() -> None:
    if not PLOT_PATH.exists():
        raise FileNotFoundError(
            f"Missing {PLOT_PATH}; run calculations/plot_dirac_potentials.py"
        )


def build_story():
    story = []

    cover_title = ParagraphStyle(
        "DiracCoverTitle",
        fontName=FONT_BOLD,
        fontSize=27,
        leading=31,
        textColor=NAVY,
        alignment=TA_LEFT,
    )
    cover_sub = ParagraphStyle(
        "DiracCoverSub",
        fontName=FONT,
        fontSize=13,
        leading=18,
        textColor=BLUE,
        alignment=TA_LEFT,
    )

    story.append(Spacer(1, 20 * mm))
    story.append(P("NOTA DI CALCOLO VERIFICABILE", "h3"))
    story.append(Spacer(1, 8 * mm))
    story.append(
        Paragraph(
            "Dirac WKB<br/>idrodinamica",
            cover_title,
        )
    )
    story.append(
        Paragraph(
            "Schwarzschild massless, gerarchia Madelung e prima estensione matriciale a Kerr",
            cover_sub,
        )
    )
    story.append(Spacer(1, 9 * mm))
    story.append(HRFlowable(width="32%", thickness=2.2, color=TEAL, hAlign="LEFT"))
    story.append(Spacer(1, 10 * mm))
    story.append(
        P(
            "<b>Risultato centrale.</b> Nel conteggio eikonale naturale di Dirac, "
            "il termine di spin connection entra a ordine 1/|kappa|, mentre il "
            "potenziale scalare di Madelung comincia a ordine 1/|kappa|^2. "
            "La gerarchia efficace non e quindi una sequenza di soli potenziali "
            "di Madelung.",
            "callout",
        )
    )
    story.append(Spacer(1, 9 * mm))
    meta = [
        ["Campo di prova", "Dirac massless su Schwarzschild"],
        ["Estensione", "Sistema radiale 2 x 2 e ricorsione proiettiva su Kerr"],
        ["Convenzioni", "G = c = hbar = 1; exp(-i omega t); Im omega < 0"],
        ["Data", "26 agosto 2026"],
    ]
    story.append(table(meta, [37 * mm, 122 * mm], header=False, font_size=8.3, leading=11))
    story.append(Spacer(1, 14 * mm))
    story.append(P("Esito operativo", "h3"))
    story.append(
        Paragraph(
            "Schwarzschild e gia abbastanza ricco per distinguere chiaramente "
            "geometria, spin, pressione quantistica e flusso QNM. Kerr non e "
            "un esercizio cosmetico: introduce frame dragging, autovalore "
            "angolare dipendente da a omega e una gerarchia naturalmente "
            "matriciale o proiettiva.",
            ParagraphStyle(
                "CoverOutcome",
                fontName=FONT_BOLD,
                fontSize=11.4,
                leading=16,
                textColor=NAVY,
            ),
        )
    )
    story.append(PageBreak())

    story += section("1. Mappa dei risultati")
    story.append(
        VerdictCard(
            159 * mm,
            "Dimostrato",
            "Partner di Schwarzschild",
            "V_+ e V_- sono ottenuti dal sistema radiale di primo ordine e verificati simbolicamente.",
            TEAL,
        )
    )
    story.append(Spacer(1, 4 * mm))
    story.append(
        VerdictCard(
            159 * mm,
            "Nuovo ordinamento",
            "Spin prima di Madelung",
            "Il simbolo subprincipale di spin e O(epsilon); Q_M comincia a O(epsilon^2).",
            BLUE,
        )
    )
    story.append(Spacer(1, 4 * mm))
    story.append(
        VerdictCard(
            159 * mm,
            "Aperto",
            "Chiusura covariante",
            "La ricorsione locale esiste; il ponte uniforme con i Lambda_n QNM e la versione bilineare covariante restano da provare.",
            GOLD,
        )
    )
    story.append(Spacer(1, 7 * mm))
    story.append(
        P(
            "La Madelung di una componente scalare disaccoppiata e utile, ma "
            "non coincide con l'idrodinamica completa del bispinore. Per "
            "quest'ultima occorrono densita, corrente e polarizzazione interna."
        )
    )
    story.append(
        table(
            [
                ["Livello", "Oggetto", "Ordine eikonale"],
                ["Principale", "h(r)^2 = f/r^2", "O(1)"],
                ["Spinoriale", "sigma s epsilon h'", "O(epsilon)"],
                ["Madelung", "-epsilon^2 A''/A", "O(epsilon^2)"],
                ["Misto", "Q_3[q_0,q_1] e successivi", "O(epsilon^3), ..."],
            ],
            [33 * mm, 85 * mm, 41 * mm],
            font_size=8,
            leading=10.5,
        )
    )

    story += section("2. Schwarzschild: dal sistema di Dirac ai partner")
    story.append(Eq("f = 1 - 2M/r,     dr*/dr = 1/f,     D* = f d/dr"))
    story.append(
        P(
            "Con K = |kappa| = j + 1/2 e una scelta dichiarata delle fasi "
            "radiali, il sistema massless e"
        )
    )
    story.append(Eq("(D* - W)F = i omega G,     (D* + W)G = i omega F"))
    story.append(Eq("W = K sqrt(f)/r"))
    story.append(P("Il disaccoppiamento fornisce esattamente"))
    story.append(Eq("F'' + (omega^2 - V_+)F = 0,     G'' + (omega^2 - V_-)G = 0"))
    story.append(Eq("V_+/- = W^2 +/- W',     W' = K sqrt(f)(3M-r)/r^3"))
    story.append(
        Eq(
            "V_+/- = K^2 f/r^2 +/- K sqrt(f)(3M-r)/r^3.                 (2.1)"
        )
    )
    story.append(
        P(
            "Per M = 1, V_+ coincide con l'equazione (37) di Cho. Il cambio "
            "di segno di kappa scambia i partner. Poiche W tende a zero ai "
            "due estremi, l'intertwiner di Darboux preserva le condizioni QNM "
            "per omega non nullo; lo spettro esatto e isospettrale."
        )
    )
    story.append(
        Eq(
            "orizzonte: Z ~ exp(-i omega r*)     |     infinito: Z ~ exp(+i omega r*)"
        )
    )
    story.append(
        P(
            "Per Im omega < 0 le autofunzioni QNM crescono spazialmente agli "
            "estremi e non sono stati L2. Una densita idrodinamica globale "
            "positiva non e quindi disponibile senza una formulazione di "
            "sistema aperto."
        )
    )
    story.append(Spacer(1, 4 * mm))
    story.append(Image(str(PLOT_PATH), width=151 * mm, height=84 * mm))
    story.append(
        P(
            "Figura 1 - I due partner standard e il potenziale stampato nel "
            "preprint arXiv:2605.28887v1. Il controllo e eseguito per M = K = 1.",
            "caption",
        )
    )

    story += section("3. Conteggio eikonale e chiusura Madelung")
    story.append(
        P(
            "La scelta decisiva e inserire epsilon nel sistema di Dirac prima "
            "di quadrarlo. Poniamo K = 1/epsilon, omega = K Omega e "
            "s = sgn(kappa), con h = sqrt(f)/r."
        )
    )
    story.append(
        Eq(
            "[-i epsilon sigma_x D* + s h sigma_y] Psi = Omega Psi.       (3.1)"
        )
    )
    story.append(
        Eq(
            "epsilon^2 Z_sigma'' + [q_0 + epsilon q_1]Z_sigma = 0,<br/>"
            "q_0 = Omega^2 - h^2,     q_1 = -sigma s h'.                 (3.2)"
        )
    )
    story.append(
        P(
            "In una regione oscillatoria reale scriviamo "
            "Z = P^(-1/2) exp[(i/epsilon) integral P dx]. La sostituzione non "
            "e approssimata:"
        )
    )
    story.append(
        Eq(
            "q_0 + epsilon q_1 = P^2 + Q_M,<br/>"
            "Q_M = -epsilon^2 (P^(-1/2))'' / P^(-1/2).                  (3.3)"
        )
    )
    story.append(
        Eq(
            "Omega^2 = P^2 + h^2 + sigma s epsilon h' + Q_M.            (3.4)"
        )
    )
    story.append(P("Per P = P_0 + epsilon P_1 + epsilon^2 P_2 + ...:"))
    story.append(
        Eq(
            "P_0 = sqrt(q_0),     P_1 = q_1/(2P_0),<br/>"
            "P_2 = 5(q_0')^2/(32 q_0^(5/2)) - q_0''/(8 q_0^(3/2)) "
            "- q_1^2/(8 q_0^(3/2))."
        )
    )
    story.append(
        Eq(
            "Q_2 = q_0''/(4q_0) - 5(q_0')^2/(16q_0^2).                  (3.5)"
        )
    )
    story.append(
        Eq(
            "Q_3 = q_1''/(4q_0) - 5q_0'q_1'/(8q_0^2) "
            "- q_1q_0''/(4q_0^2) + 5q_1(q_0')^2/(8q_0^3).              (3.6)"
        )
    )
    story.append(
        P(
            "Il termine Q_3 e il primo contributo misto spin-gradiente. "
            "La presenza del simbolo subprincipale q_1 rompe la parita della "
            "serie di Madelung completa, anche se il funzionale porta un "
            "prefattore epsilon^2."
        )
    )
    story.append(
        Eq(
            "2P_0 P_n = q_n + [epsilon^(n-2)]F[P] "
            "- sum_(j=1)^(n-1) P_j P_(n-j),<br/>"
            "F[P] = 3/4(P'/P)^2 - 1/2 P''/P.                           (3.7)"
        )
    )
    story.append(
        P(
            "La fase di spin del primo ordine si integra localmente: "
            "integral P_1 dr* = -(sigma s/2) arcsin(h/Omega), per un ramo "
            "fissato e q_0 > 0."
        )
    )
    story.append(
        KeepTogether([P(
            "<b>Limite.</b> Q_2, Q_3 e i successivi contengono potenze inverse "
            "di q_0 e divergono ai turning point. Non possono essere inseriti "
            "direttamente nel massimo della barriera QNM: serve la forma "
            "normale uniforme di Iyer-Will.",
            "callout",
        )])
    )

    story += section("4. Idrodinamica del bispinore")
    story.append(
        P(
            "La corrente fisica emerge dalla Hamiltoniana 2 x 2, non da una "
            "sola componente. Nella base di (3.1):"
        )
    )
    story.append(
        Eq(
            "rho = Psi^dagger Psi,     j = Psi^dagger sigma_x Psi,<br/>"
            "partial_t rho + partial_* j = 0.                           (4.1)"
        )
    )
    story.append(
        Eq(
            "(p sigma_x + s h sigma_y)chi = Omega chi,<br/>"
            "Omega^2 = p^2 + h^2,     j = rho p/Omega.                  (4.2)"
        )
    )
    story.append(
        P(
            "Dopo una rotazione unitaria costante, H = -i sigma_3 D* + W "
            "sigma_1. Le quattro variabili reali di Bloch chiudono il sistema:"
        )
    )
    story.append(
        Eq(
            "N = Psi^dagger Psi,  J = Psi^dagger sigma_3 Psi,<br/>"
            "Sigma_r = Psi^dagger sigma_1 Psi,  Pi_r = Psi^dagger sigma_2 Psi"
        )
    )
    story.append(Eq("N^2 = J^2 + Sigma_r^2 + Pi_r^2."))
    story.append(
        Eq(
            "N' = 2W Pi_r - 2 omega_I J,     J' = -2 omega_I N,<br/>"
            "Sigma_r' = 2 omega_R Pi_r,      Pi_r' = 2WN - 2 omega_R Sigma_r.  (4.3)"
        )
    )
    story.append(
        P(
            "Per omega reale, J e costante. Ponendo N = J cosh(eta), "
            "Sigma_r = J sinh(eta)cos(beta), Pi_r = J sinh(eta)sin(beta):"
        )
    )
    story.append(
        Eq(
            "eta' = 2W sin(beta),     beta' = -2omega + 2W coth(eta)cos(beta)."
        )
    )
    story.append(
        P(
            "Questo e il candidato idrodinamico radiale piu fedele: N e J "
            "sono densita e flusso, eta e beta conservano coerenza e "
            "polarizzazione. Le densita dei partner sono "
            "|Z_sigma|^2 = (N + sigma Pi_r)/2."
        )
    )
    story.append(
        P(
            "Per una QNM, J' = -2 Im(omega) N. Il modo completo, incluso il "
            "fattore temporale, soddisfa ancora la continuita; il profilo "
            "stazionario isolato descrive invece un flusso aperto."
        )
    )

    story += section("5. Photon sphere e benchmark WKB3")
    story.append(
        Eq(
            "r_peak = 3M - sigma s sqrt(3) M/(2K) + O(K^-2),<br/>"
            "V_peak = K^2/(27M^2) + 1/(108M^2) + O(K^-1).              (5.1)"
        )
    )
    story.append(
        Eq(
            "M omega_(K,n) = [K - i(n+1/2)]/[3 sqrt(3)] + O(K^-1).      (5.2)"
        )
    )
    story.append(
        P(
            "Il coefficiente completo O(K^-1) richiede le correzioni uniformi "
            "e non segue dal solo valore del massimo. Il controllo numerico "
            "WKB di terzo ordine, per M = 1 e n = 0, e:"
        )
    )
    story.append(
        table(
            [
                ["K", "r_peak/M", "M omega WKB3"],
                ["1", "2.420483", "0.176452 - 0.100109 i"],
                ["2", "2.617914", "0.378627 - 0.096542 i"],
                ["3", "2.727147", "0.573685 - 0.096324 i"],
                ["4", "2.790186", "0.767194 - 0.096276 i"],
                ["5", "2.830203", "0.960215 - 0.096256 i"],
            ],
            [24 * mm, 47 * mm, 88 * mm],
            font_size=8.2,
            leading=10.8,
        )
    )
    story.append(
        P(
            "I valori riproducono la tabella di Cho. Per K = 1, il continued "
            "fraction di Jing fornisce M omega = 0.182963 - 0.0969825 i, "
            "mostrando il limite quantitativo della WKB3 al modo angolare piu basso."
        )
    )

    story += section("6. Kerr: gerarchia matriciale e proiettiva")
    story.append(
        Eq(
            "Delta = r^2 - 2Mr + a^2,     dr*/dr = (r^2+a^2)/Delta,<br/>"
            "v = am/(r^2+a^2),     W = lambda sqrt(Delta)/(r^2+a^2)."
        )
    )
    story.append(
        Eq(
            "(D* + i[omega-v])P_- = W P_+,<br/>"
            "(D* - i[omega-v])P_+ = W P_-.                              (6.1)"
        )
    )
    story.append(
        Eq(
            "[v I + i sigma_z D* + W sigma_y]P = omega P.               (6.2)"
        )
    )
    story.append(
        P(
            "Con m = K mu, lambda = K Lambda(a Omega,mu), omega = K "
            "Omega e epsilon = 1/K, il simbolo principale produce"
        )
    )
    story.append(
        Eq(
            "(Omega-v_0)^2 = p^2 + w_0^2,<br/>"
            "j = -P^dagger sigma_z P = rho p/(Omega-v_0).               (6.3)"
        )
    )
    story.append(
        P(
            "v_0 e il frame dragging. All'orizzonte Omega-v_0 diventa la "
            "frequenza co-rotante omega-m Omega_H; all'infinito v_0 e w_0 "
            "tendono a zero."
        )
    )
    story.append(
        P(
            "Per z = P_+/P_- e kappa_r = Omega - a mu/(r^2+a^2), "
            "il sistema ammette una Riccati proiettiva esatta:"
        )
    )
    story.append(
        Eq(
            "epsilon z' = w_0 + 2i kappa_r z - w_0 z^2.                 (6.4)"
        )
    )
    story.append(
        Eq(
            "k = sqrt(kappa_r^2-w_0^2),<br/>"
            "z_0 = i(kappa_r + tau k)/w_0,     tau = +/-1,"
        )
    )
    story.append(
        Eq(
            "z_n = [z_(n-1)' + w_0 sum_(j=1)^(n-1) z_j z_(n-j)]/"
            "[-2i tau k].                                               (6.5)"
        )
    )
    story.append(
        Eq(
            "(log P_-)' = i tau k/epsilon + w_0 z_1 + epsilon w_0 z_2 + ..."
        )
    )
    story.append(
        P(
            "Questa e gia una gerarchia all-order che conserva la "
            "polarizzazione. Un Q_eff scalare esiste solo dopo una scelta di "
            "componente, gauge e ramo; non e in generale canonico."
        )
    )
    story += subsection("Perche il superpotenziale scalare di Kerr e delicato")
    story.append(
        Eq(
            "d xhat/dr = K(r)/(omega Delta),<br/>"
            "U = lambda omega sqrt(Delta)/K(r),     V_+/- = U^2 +/- dU/dxhat."
        )
    )
    story += bullets(
        [
            "xhat e U dipendono dall'autovalore complesso omega.",
            "lambda dipende dal problema angolare e da a omega.",
            "K(r) = 0 genera una singolarita della trasformazione.",
            "La WKB scalare di Kerr e gia disponibile fino al sesto ordine; la novita plausibile e la gerarchia matriciale/proiettiva.",
            "Per Dirac neutro su Kerr non si ha superradianza classica.",
        ]
    )

    story += section("7. Letteratura e controllo del precedente diretto")
    story.append(
        P(
            "Takabayasi fornisce la cornice bilineare originaria. Matos, "
            "Gallegos e Chavanis (2022) hanno gia pubblicato una "
            "rappresentazione idrodinamica di Dirac e Weyl in spazio-tempo "
            "curvo: l'idea generale non e quindi nuova. Il contributo "
            "distintivo deve riguardare l'ordinamento WKB, la polarizzazione, "
            "i QNM e Kerr."
        )
    )
    story.append(
        P(
            "Il preprint Meza-Dominguez-Matos, arXiv:2605.28887v1 "
            "(27 maggio 2026), e direttamente sovrapposto: applica una "
            "formulazione chiral-idrodinamica a Schwarzschild, QNM e greybody "
            "factors. Deve essere citato, ma la v1 richiede controlli indipendenti.",
            "callout",
        )
    )
    audit = [
        ["Controllo", "Esito nella v1"],
        ["Limite massless", "Le velocita (3.10)-(3.11) contengono 1/m."],
        [
            "Correnti chirali",
            "La (4.20) le conserva separatamente anche per m non nullo; il termine di massa trasferisce chiralita.",
        ],
        [
            "Potenziale (7.36)",
            "Non coincide con W^2 +/- W'; per M = kappa = 1 ha r_max = 3.85313 invece di 2.42048.",
        ],
        [
            "Condizione QNM",
            "La (10.41) contiene '?'; la (10.43) divide per sqrt[f(r_s)] = 0.",
        ],
        [
            "Ramo outgoing",
            "exp(+i omega r*) cresce, non decade, per Im omega < 0.",
        ],
        [
            "Greybody (11.47)",
            "I due membri stampati non sono algebricamente uguali.",
        ],
    ]
    story.append(
        KeepTogether(
            [table(audit, [43 * mm, 116 * mm], font_size=7.25, leading=9.4)]
        )
    )
    story.append(
        P(
            "Questi rilievi riguardano le equazioni stampate nella versione "
            "v1 e non costituiscono un giudizio definitivo sul programma. "
            "La baseline tecnica resta il sistema standard di Chandrasekhar "
            "e i partner verificati di Cho/Jing."
        )
    )

    story += section("8. Stato logico e prossimi calcoli")
    status = [
        ["Affermazione", "Stato"],
        ["Partner Schwarzschild e corrente radiale", "Dimostrati"],
        ["Gerarchia P_n, Q_2, Q_3", "Identita locale; non uniforme ai turning point"],
        ["Sistema reale di Bloch", "Esatto nel canale radiale"],
        ["Ricorsione proiettiva Kerr", "Formale all-order, verificata simbolicamente ai primi ordini"],
        ["Mappa con Lambda_n QNM", "Aperta; richiede uniformizzazione"],
        ["Fluido covariante completo", "Aperto; servono bilineari e trasporto di spin"],
    ]
    story.append(table(status, [78 * mm, 81 * mm], font_size=7.7, leading=10.2))
    story += subsection("Sequenza raccomandata")
    story.append(
        table(
            [
                ["1", "Uniformizzare Schwarzschild al massimo mantenendo q_1 e verificare Lambda_2, Lambda_3."],
                ["2", "Ricostruire J^mu dalle soluzioni WKB F,G fino a O(epsilon^2)."],
                ["3", "Testare l'isospectralita dei partner dopo troncamento e Pade."],
                ["4", "Su Kerr, calcolare autofibre, connessione di Berry e dipendenza angolare."],
                ["5", "Confrontare K = 1,2 con continued fractions prima dell'interpretazione fisica."],
            ],
            [10 * mm, 149 * mm],
            header=False,
            font_size=7.55,
            leading=9.7,
        )
    )

    story += section("Riferimenti primari essenziali")
    refs = [
        (
            1,
            "Chandrasekhar (1976), separazione dell'equazione di Dirac in Kerr.",
            "https://doi.org/10.1098/rspa.1976.0090",
        ),
        (
            2,
            "Page (1976), Dirac in Kerr-Newman.",
            "https://doi.org/10.1103/PhysRevD.14.1509",
        ),
        (
            3,
            "Cho (2003), QNM di Dirac massless e massive su Schwarzschild.",
            "https://arxiv.org/abs/gr-qc/0303078",
        ),
        (
            4,
            "Jing (2005), continued fractions per QNM di Dirac su Schwarzschild.",
            "https://arxiv.org/abs/gr-qc/0502023",
        ),
        (
            5,
            "Schutz-Will (1985), WKB di barriera per QNM.",
            "https://doi.org/10.1086/184453",
        ),
        (
            6,
            "Iyer-Will (1987), WKB di ordine superiore.",
            "https://doi.org/10.1103/PhysRevD.35.3621",
        ),
        (
            7,
            "Takabayasi (1957), idrodinamica relativistica della materia di Dirac.",
            "https://doi.org/10.1143/PTPS.4.2",
        ),
        (
            8,
            "Matos-Gallegos-Chavanis (2022), idrodinamica di Dirac/Weyl in spazio-tempo curvo.",
            "https://doi.org/10.1140/epjc/s10052-022-10853-5",
        ),
        (
            9,
            "Oancea-Kumar (2023), analisi semiclassica covariante dei campi di Dirac.",
            "https://doi.org/10.1103/PhysRevD.107.044029",
        ),
        (
            10,
            "Carlson-Cornell-Jordan (2012), WKB fermionica su Kerr.",
            "https://arxiv.org/abs/1201.3267",
        ),
        (
            11,
            "Meza-Dominguez-Matos (2026), precedente chiral-idrodinamico diretto.",
            "https://arxiv.org/abs/2605.28887",
        ),
        (
            12,
            "Unruh (1973), neutrini in Kerr e assenza di superradianza fermionica classica.",
            "https://doi.org/10.1103/PhysRevLett.31.1265",
        ),
    ]
    for n, citation, url in refs:
        story.append(P(f"[{n}] {citation} {link(url, 'link')}", "ref"))

    story += section("Appendice: file riproducibili")
    files = [
        ["File", "Contenuto verificato"],
        [
            "calculations/dirac_schwarzschild_wkb.py",
            "Partner, P_0/P_1/P_2, Q_2/Q_3, massimo eikonale.",
        ],
        [
            "calculations/dirac_schwarzschild_wkb3.py",
            "Massimo esatto e frequenze WKB3 di Cho.",
        ],
        [
            "calculations/dirac_kerr_projective_wkb.py",
            "Riccati Kerr e ricorsione z_0/z_1/z_2.",
        ],
        [
            "calculations/dirac_wkb_hydrodynamics_notes.md",
            "Nota estesa con formule, limiti e audit bibliografico.",
        ],
    ]
    story.append(table(files, [84 * mm, 75 * mm], font_size=7.15, leading=9.6))
    story.append(Spacer(1, 7 * mm))
    story.append(HRFlowable(width="100%", thickness=1.2, color=TEAL))
    story.append(
        P(
            "Conclusione: su Schwarzschild la prima correzione e spinoriale, "
            "non Madelung; su Kerr la struttura naturale e matriciale o "
            "proiettiva. Questa formulazione resta valida sia che emerga, sia "
            "che fallisca, un Q_eff scalare canonico.",
            "callout",
        )
    )
    return story


def main() -> None:
    make_potential_plot()
    doc = DiracNoteDocTemplate(str(PDF_PATH))
    doc.build(build_story())
    print(PDF_PATH)


if __name__ == "__main__":
    main()
