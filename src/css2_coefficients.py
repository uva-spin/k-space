"""Massless SU(3) quark CSS2 ingredients, a=alpha_s/(4*pi).

Conventions: beta=d a/d ln(mu)=-2 sum beta_n a^(n+2);
Gamma=gamma_K/2; gamma=gamma_j; D=-K, d=d(a, L_b=0).
Source/equation locators and external-file hashes: manifests/css2_sources.json.
This is N3LL evolution with strict NLO boundaries, NOT an N3LL-prime payload.
"""
from math import pi
from scipy.special import zeta
from resummed_operator import Coefficients

CF, CA, TF = 4/3, 3., .5


def qcd_coefficients(nf=5):
    if isinstance(nf, bool) or nf not in (3, 4, 5, 6):
        raise ValueError('Supported fixed massless flavor counts are 3,4,5,6')
    z2,z3,z4,z5 = (float(zeta(n)) for n in (2,3,4,5))
    # van Ritbergen, Vermaseren, Larin hep-ph/9701390v1, SU(3) result.
    beta = (11-2*nf/3, 102-38*nf/3,
            2857/2-5033*nf/18+325*nf**2/54,
            149753/6+3564*z3-(1078361/162+6508*z3/27)*nf
            +(50065/162+6472*z3/81)*nf**2+1093*nf**3/729)
    # Collins--Rogers (59), gamma_K/2.
    cusp = [4*CF, CF*(CA*(268/9-4*pi**2/3)-40*nf/9),
            CF*(-16*nf**2/27+CA*nf*(-836/27+80*pi**2/27-112*z3/3)
                +CA**2*(490/3-536*pi**2/27+44*pi**4/45+88*z3/3)
                +CF*nf*(-110/3+32*z3))]
    # von Manteuffel--Panzer--Schabinger 2002.04617v2, eq:QCDcusp.
    # Quartic invariants are normalized by N_F=3, NOT N_A=8.
    dFA_NF, dFF_NF = 5/2, 5/36
    cusp.append(
        nf**3*CF*(64*z3/27-32/81)
        +nf**2*CA*CF*(-224*z2**2/15+2240*z3/27-608*z2/81+923/81)
        +nf**2*CF**2*(64*z2**2/5-640*z3/9+2392/81)
        +nf*CA**2*CF*(2096*z5/9+448*z3*z2/3-352*z2**2/15
                        -23104*z3/27+20320*z2/81-24137/81)
        +nf*CA*CF**2*(160*z5-128*z3*z2-352*z2**2/5+3712*z3/9+440*z2/3-34066/81)
        +nf*CF**3*(-320*z5+592*z3/3+572/9)
        +nf*dFF_NF*(-1280*z5/3-256*z3/3+256*z2)
        +dFA_NF*(-384*z3**2-7936*z2**3/35+3520*z5/3+128*z3/3-128*z2)
        +CA**3*CF*(-16*z3**2-20032*z2**3/105-3608*z5/9-352*z3*z2/3
                      +3608*z2**2/5+20944*z3/27-88400*z2/81+84278/81))
    # Collins--Rogers (58).
    gamma = (6*CF,
        CF**2*(3-4*pi**2+48*z3)+CF*CA*(961/27+11*pi**2/3-52*z3)
            +CF*nf*(-130/27-2*pi**2/3),
        CF**2*nf*(-2953/27+26*pi**2/9+28*pi**4/27-512*z3/9)
        +CF*nf**2*(-4834/729+20*pi**2/27+16*z3/27)
        +CF**3*(29+6*pi**2+16*pi**4/5+136*z3-32*pi**2*z3/3-480*z5)
        +CA**2*CF*(139345/1458+7163*pi**2/243+83*pi**4/45-7052*z3/9+88*pi**2*z3/9+272*z5)
        +CA*CF*nf*(17318/729-2594*pi**2/243-22*pi**4/45+1928*z3/27)
        +CA*CF**2*(151/2-410*pi**2/9-494*pi**4/135+1688*z3/3+16*pi**2*z3/3+240*z5))
    # Li--Zhu 1604.01404v4 (9); CR VIII.3 and Appendix B: K=gamma_r.
    gr1=CF*CA*(28*z3-808/27)+112*CF*nf/27
    gr2=(CF*CA**2*(-176*z3*z2/3+6392*z2/81+12328*z3/27+154*z4/3-192*z5-297029/729)
         +CF*CA*nf*(-824*z2/81-904*z3/27+20*z4/3+62626/729)
         +CF*nf**2*(-32*z3/9-1856/729)
         +CF**2*nf*(-304*z3/9-16*z4+1711/27))
    return Coefficients(
        tuple(beta), tuple(cusp), tuple(gamma), (0., -gr1, -gr2),
        description=f'Physical massless SU(3), fixed nf={nf}; N3LL evolution, strict NLO boundaries',
        physical_n3llprime_payload=False,
        scheme='CSS2 D=-K; a=alpha_s/(4*pi); derivatives ln(mu)',
        evolution_accuracy='N3LL', hard_and_ope_order='strict NLO', nf=nf,
        source_manifest='manifests/css2_sources.json')
