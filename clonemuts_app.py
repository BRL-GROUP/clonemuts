import streamlit as st
import math as m

# ── stringsandbiophysics core ──────────────────────────────────────────────

R = 1.9865

deltaH = {
    "AA":-7.9,"TT":-7.9,"AT":-7.2,"TA":-7.2,"AC":-8.4,"GT":-8.4,
    "AG":-7.8,"CT":-7.8,"TC":-8.2,"GA":-8.2,"TG":-8.5,"CA":-8.5,
    "CC":-8.0,"GG":-8.0,"CG":-10.6,"GC":-9.8,"initA/T":2.3,"initG/C":0.1
}
deltaS = {
    "AA":-22.2,"TT":-22.2,"AT":-20.4,"TA":-21.3,"AC":-22.4,"GT":-22.4,
    "AG":-21.0,"CT":-21.0,"TC":-22.2,"GA":-22.2,"TG":-22.7,"CA":-22.7,
    "CC":-19.9,"GG":-19.9,"CG":-27.2,"GC":-24.4,"initA/T":4.1,"initG/C":-2.8
}

workingConcs = {"oligoconc":500.0,"Mgconc":2.0,"MonoValentconc":50.0,"dNTPsconc":0.8}
oligoconc      = workingConcs['oligoconc']
Mgconc         = workingConcs['Mgconc']
MonoValentconc = workingConcs['MonoValentconc']
dNTPsconc      = workingConcs['dNTPsconc']

MV     = MonoValentconc*10**-3.0
LogMV  = m.log(MV)
Log2MV = LogMV*LogMV
Log3MV = LogMV*LogMV*LogMV

if Mgconc > 0.0:
    MgKa   = 3.0*10**4.0
    FreeMg = (-(MgKa*(dNTPsconc-Mgconc)*10**-3+1)+m.sqrt((MgKa*(dNTPsconc-Mgconc)*10**-3+1)**2+4*MgKa*Mgconc*10**-3))/(2*MgKa)
    LogMg  = m.log(FreeMg)
    Log2Mg = LogMg*LogMg
else:
    FreeMg = 0.0

Ratio = m.sqrt(FreeMg)/(MonoValentconc*10**-3.0)

def melting(melt):
    oligolength = len(melt)
    if oligolength < 8:
        return -100.0
    Cc = melt.count('C'); Gc = melt.count('G')
    FracGC = float(Cc+Gc)/float(oligolength)
    dH = dS = 0.0
    for key in deltaH:
        pos = keycount = 0
        while True:
            try:
                idx = melt.index(key, pos); keycount += 1; pos = idx+1
            except ValueError:
                break
        dH += keycount*deltaH[key]; dS += keycount*deltaS[key]
    for end in (melt[0], melt[-1]):
        tag = "initA/T" if end in "AT" else "initG/C"
        dH += deltaH[tag]; dS += deltaS[tag]
    TmBase = (dH*1000.0)/(dS + R*m.log(oligoconc*10.0**-9.0))
    if Ratio < 0.22:
        inVTm = (1.0/TmBase)+((4.29*FracGC-3.95)*LogMV+0.940*Log2MV)*10.0**-5
    elif Ratio <= 6.0:
        af = 3.92*(0.843-0.352*m.sqrt(MV)*LogMV)
        df = 1.42*(1.279-10**-3*(4.03*LogMV+8.03*Log2MV))
        gf = 8.31*(0.486-0.258*LogMV+5.25*10**-3*Log3MV)
        inVTm = (1.0/TmBase+(af-0.911*LogMg+FracGC*(6.26+df*LogMg)+1.0/(2.0*(oligolength-1.0))*(-48.2+52.5*LogMg+gf*Log2Mg))*10**-5)
    else:
        inVTm = (1.0/TmBase+(3.92-0.911*LogMg+FracGC*(6.26+1.42*LogMg)+1.0/(2.0*(oligolength-1.0))*(-48.2+52.5*LogMg+8.31*Log2Mg))*10**-5)
    return 1.0/inVTm - 273.15

_table = str.maketrans('ACBDGHK\nMNSRUTWVYacbdghkmnsrutwvy',
                        'TGVHCDM\nKNSYAAWBRTGVHCDMKNSYAAWBR')

def rvscomp(seq):
    return (''.join(seq)).translate(_table)[::-1]

def seqbreakfwd(dna, Tg):
    melt = dna[0:8]; Tm = melting(melt)
    for i in range(1, len(dna)-8):
        if Tm >= Tg: return melt, Tm
        melt = dna[0:8+i]; Tm = melting(melt)
    return melt, Tm

def seqbreakrvs(dna, Tg):
    melt = dna[len(dna)-8:]; Tm = melting(melt)
    for i in range(1, len(dna)-8):
        if Tm >= Tg: return melt, Tm
        melt = dna[len(dna)-8-i:]; Tm = melting(melt)
    return melt, Tm

# ── Cloneit logic ──────────────────────────────────────────────────────────

def cloneit(insert, upstream, downstream, Tg):
    ifwd,  Tmifwd  = seqbreakfwd(insert, Tg)
    vrvsovcomp, Tm5ov = seqbreakrvs(upstream, Tg)
    vrvsov = rvscomp(vrvsovcomp)
    split  = upstream[0:len(upstream)-len(vrvsovcomp)]
    vrvscomp2, Tmvrvs = seqbreakrvs(split, Tg)
    vrvs   = rvscomp(vrvscomp2)
    ovifwd = vrvsovcomp + ifwd
    ovvrvs = vrvsov + vrvs
    irvscomp, Tmirvs = seqbreakrvs(insert, Tg)
    irvs   = rvscomp(irvscomp)
    vfwdov, Tm3ov = seqbreakfwd(downstream, Tg)
    vfwdovcomp = rvscomp(vfwdov)
    split2 = downstream[len(vfwdov):]
    vfwd, Tmvfwd = seqbreakfwd(split2, Tg)
    ovvfwd = vfwdov + vfwd
    ovirvs = vfwdovcomp + irvs

    seqs = [ovvfwd, ovvrvs, ovifwd, ovirvs]
    maxlen = max(len(s) for s in seqs)
    fs = "%-9s %-" + str(maxlen) + "s  Length: %3s  Tm: %s\n"
    fl = "%-9s %-" + str(maxlen) + "s  Length: %3s\n"

    out  = fs % ("vecfwd:",   vfwd,    len(vfwd),    f"{Tmvfwd:.1f} C")
    out += fl % ("ovvecfwd:", ovvfwd,  len(ovvfwd))
    out += fs % ("vecrvs:",   vrvs,    len(vrvs),    f"{Tmvrvs:.1f} C")
    out += fl % ("ovvecrvs:", ovvrvs,  len(ovvrvs))
    out += fs % ("insfwd:",   ifwd,    len(ifwd),    f"{Tmifwd:.1f} C")
    out += fl % ("ovinsfwd:", ovifwd,  len(ovifwd))
    out += fs % ("insrvs:",   irvs,    len(irvs),    f"{Tmirvs:.1f} C")
    out += fl % ("ovinsrvs:", ovirvs,  len(ovirvs))
    out += f"3' overlap Tm:  {Tm3ov:.1f} C\n"
    out += f"5' overlap Tm:  {Tm5ov:.1f} C"
    return out

# ── Mutit logic ────────────────────────────────────────────────────────────

def mutit(mutation, upstream, downstream, Tg):
    melt = mutation; vrvs = upstream; vfwd = downstream
    Tm   = melting(melt)
    for i in range(1, len(upstream)+len(downstream)):
        if Tm >= Tg: break
        melt = upstream[len(vrvs)-1:len(vrvs)] + melt
        vrvs = upstream[0:len(vrvs)-1]; Tm = melting(melt)
        if Tm >= Tg: break
        melt = melt + downstream[len(downstream)-len(vfwd):len(downstream)-len(vfwd)+1]
        vfwd = downstream[len(downstream)-len(vfwd)+1:]; Tm = melting(melt)
    Tmov     = Tm
    ovfwd    = melt
    vfwd_seq, Tmvfwd = seqbreakfwd(vfwd, Tg)
    vrvscomp, Tmvrvs = seqbreakrvs(vrvs, Tg)
    vrvs_seq         = rvscomp(vrvscomp)
    ovvfwd           = ovfwd + vfwd_seq
    ovrvs            = rvscomp(ovfwd)
    ovvrvs           = ovrvs + vrvs_seq

    maxlen = max(len(ovvfwd), len(ovvrvs))
    fs = "%-9s %-" + str(maxlen) + "s  Length: %3s  Tm: %s\n"
    fl = "%-9s %-" + str(maxlen) + "s  Length: %3s\n"

    out  = fs % ("vecfwd:",   vfwd_seq, len(vfwd_seq), f"{Tmvfwd:.1f} C")
    out += fl % ("ovvecfwd:", ovvfwd,   len(ovvfwd))
    out += fs % ("vecrvs:",   vrvs_seq, len(vrvs_seq), f"{Tmvrvs:.1f} C")
    out += fl % ("ovvecrvs:", ovvrvs,   len(ovvrvs))
    out += f"Overlap Tm:  {Tmov:.1f} C"
    return out

# ── Helpers ────────────────────────────────────────────────────────────────

NUCLEOTIDES = set('ACBDGHKMNSUTWVY')

def clean(seq):
    return seq.replace(" ","").replace("\n","").replace("\r","").upper()

def valid(seq):
    return all(c in NUCLEOTIDES for c in seq)

# ── Streamlit UI ───────────────────────────────────────────────────────────

st.set_page_config(page_title="CloneMuts", page_icon="🧬", layout="wide")
st.title("🧬 CloneMuts")
st.caption("Primer design for cloning and site-directed mutagenesis")

tab_clone, tab_mut = st.tabs(["Cloning", "Mutagen"])

# ── Cloning tab ────────────────────────────────────────────────────────────
with tab_clone:
    st.subheader("Cloning Primer Design")
    col1, col2, col3 = st.columns(3)
    with col1:
        insert_c     = st.text_area("Insert sequence", height=120, key="insert_c")
    with col2:
        upstream_c   = st.text_area("Upstream vector sequence (5')", height=120, key="up_c")
    with col3:
        downstream_c = st.text_area("Downstream vector sequence (3')", height=120, key="dn_c")

    tm_c = st.number_input("Target Tm (°C)", min_value=30.0, max_value=80.0,
                            value=52.0, step=0.5, key="tm_c")

    if st.button("Design Cloning Primers", type="primary", key="btn_clone"):
        ic = clean(insert_c); uc = clean(upstream_c); dc = clean(downstream_c)
        errors = []
        if len(ic) < 8:  errors.append("Insert sequence too short (min 8 bp)")
        if len(uc) < 8:  errors.append("Upstream sequence too short (min 8 bp)")
        if len(dc) < 8:  errors.append("Downstream sequence too short (min 8 bp)")
        if ic and not valid(ic): errors.append("Insert contains non-nucleotide characters")
        if uc and not valid(uc): errors.append("Upstream contains non-nucleotide characters")
        if dc and not valid(dc): errors.append("Downstream contains non-nucleotide characters")

        if errors:
            for e in errors: st.error(e)
        else:
            try:
                out = cloneit(ic, uc, dc, float(tm_c))
                st.success("Primers designed successfully!")
                st.text_area("Results — select all and copy (⌘A / Ctrl+A)",
                             value=out, height=300, key="out_clone")
            except Exception as e:
                st.error(f"Could not design primers: {e}. Check sequences are long enough.")

# ── Mutagen tab ────────────────────────────────────────────────────────────
with tab_mut:
    st.subheader("Mutagenesis Primer Design")
    col1, col2, col3 = st.columns(3)
    with col1:
        mutation_m   = st.text_area("Mutation sequence (leave empty to excise)", height=120, key="mut_m")
    with col2:
        upstream_m   = st.text_area("Upstream vector sequence (5')", height=120, key="up_m")
    with col3:
        downstream_m = st.text_area("Downstream vector sequence (3')", height=120, key="dn_m")

    tm_m = st.number_input("Target Tm (°C)", min_value=30.0, max_value=80.0,
                            value=52.0, step=0.5, key="tm_m")

    mc_preview = clean(mutation_m)
    if 40 < len(mc_preview) <= 60: st.warning("Long mutation sequence — primers may be unwieldy")
    if len(mc_preview) > 60:       st.error("Mutation sequence too long (max 60 bp)")

    if st.button("Design Mutagenesis Primers", type="primary", key="btn_mut"):
        mc = clean(mutation_m); um = clean(upstream_m); dm = clean(downstream_m)
        errors = []
        if len(um) < 8:  errors.append("Upstream sequence too short (min 8 bp)")
        if len(dm) < 8:  errors.append("Downstream sequence too short (min 8 bp)")
        if len(mc) > 60: errors.append("Mutation sequence too long (max 60 bp)")
        if mc and not valid(mc): errors.append("Mutation contains non-nucleotide characters")
        if um and not valid(um): errors.append("Upstream contains non-nucleotide characters")
        if dm and not valid(dm): errors.append("Downstream contains non-nucleotide characters")

        if errors:
            for e in errors: st.error(e)
        else:
            try:
                out = mutit(mc, um, dm, float(tm_m))
                st.success("Primers designed successfully!")
                st.text_area("Results — select all and copy (⌘A / Ctrl+A)",
                             value=out, height=200, key="out_mut")
            except Exception as e:
                st.error(f"Could not design primers: {e}. Check sequences are long enough.")
