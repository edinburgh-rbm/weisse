 /* 
  *  MOIS Examples: Mammalian Circadian Clock (Forger-Peskin 2002)
  *  Copyright (C) 2014 University of Edinburgh School of Informatics
  * 
  *  This program is free software: you can redistribute it and/or modify
  *  it under the terms of the GNU General Public License as published by
  *  the Free Software Foundation, either version 3 of the License, or
  *  (at your option) any later version.
  * 
  *  This program is distributed in the hope that it will be useful,
  *  but WITHOUT ANY WARRANTY; without even the implied warranty of
  *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  *  GNU General Public License for more details.
  * 
  *  You should have received a copy of the GNU General Public License
  *  along with this program. If not, see <http://www.gnu.org/licenses/>.
  */
package uk.ac.ed.inf.mois.clocks

import uk.ac.ed.inf.mois.{Model, ODE, ProcessGroup}
import uk.ac.ed.inf.mois.sched.NaiveScheduler
import spire.implicits._
import uk.ac.ed.inf.mois.implicits._

class MammalianCircadianClock(L: Double)
  extends ODE {

 /* Global variables */

 val tmc = Double("c:tmc") default(0.42)
 tmc annotate("description", "Preparation and nuclear export of all mRNA")
 tmc annotate("units", "1/h")

 val umR = Double("c:umR") default(0.3)
 umR annotate("description", "Degradation of CRY1 and CRY2 mRNA")
 umR annotate("units", "1/h")

 val up = Double("c:up") default(3.39)
 up annotate("description", "Degradation of PERp unbound to CRY")
 up annotate("units", "1/h")

 val ar = Double("c:ar") default(0.27)
 ar annotate("description", "Binding of PER1 and PER2 to CRY1 and CRY2")
 ar annotate("units", "1/nM 1/h")

 val dr = Double("c:dr") default(0.35)
 dr annotate("description", "Unbinding of PER1 and PER2 from CRY1 and CRY2")
 dr annotate("units", "1/h")

 val nl = Double("c:nl") default(2.31)
 nl annotate("description", "Nuclear localisation of PERp and bound proteins")
 nl annotate("units", "1/h")

 val ne = Double("c:ne") default(0.71)
 ne annotate("description", "Nuclear export of PER and boun proteins")
 ne annotate("units", "1/h")

 val Nf = Double("c:Nf") default(115.76)
 Nf annotate("description", "Ratio of nuclear to cytoplasmic compartment volume")


 /* Probabilities of promoter binding */

 val G = Double("c:G") default(0.884)
 G annotate("description", "Probability of CRY bound to promoter")

 val GRv = Double("c:GRv") default(0)
 GRv annotate("description", "Probability of REVERBa bound to promoter")

 val bin = Double("c:bin") default(1476.52)
 bin annotate("description", "Binding of CRY to promoter in nucleus")
 bin annotate("units", "1/nM 1/h")

 val unbin = Double("c:unbin") default(23.78)
 unbin annotate("description", "Unbinding of CRY from promoter in nucleus")
 unbin annotate("units", "1/h")

 val binRv = Double("c:binRv") default(0.13)
 binRv annotate("description", "Normalised binding of REVERBa to promoter in nucleus")
 binRv annotate("units", "1/nM 1/h")

 val unbinRv = Double("c:unbinRv") default(21.76)
 unbinRv annotate("description", "Normalised unbinding of REVERBa from promoter in nucleus")
 unbinRv annotate("units", "1/h")


 /* Kinase */

 val Cn = Double("c:Cn") default(0.31)
 Cn annotate("description", "Nuclear kinase")
 Cn annotate("units", "nM")

 val Ct = Double("c:Ct") default(7.71)
 Ct annotate("description", "Total kinase concentration")
 Ct annotate("units", "nM")

 val ac = Double("c:ac") default(0.47)
 ac annotate("description", "Binding of PER1 and PER2 to kinases")
 ac annotate("units", "1/nM 1/h")

 val dc = Double("c:dc") default(5.09)
 dc annotate("description", "Unbinding of PER1 and PER2 from kinases")
 dc annotate("units", "1/h")

 val hoo = Double("c:hoo") default(0.29)
 hoo annotate("description", "Initial phosphorylation of PER1")
 hoo annotate("units", "1/h")

 val hot = Double("c:hot") default(0.09)
 hot annotate("description", "Initial phosphorylation of PER2")
 hot annotate("units", "1/h")

 val hto = Double("c:hto") default(1.45)
 hto annotate("description", "Phosphorylation preventing nuclear localisation of PER1")
 hto annotate("units", "1/h")

 val htt = Double("c:htt") default(0)
 htt annotate("description", "Phosphorylation preventing nuclear localisation of PER2")
 ht annotate("units", "1/h")


 /* PER1 Gene */

 val MnPo = Double("c:MnPo") default(0.039)
 MnPo annotate("description", "Nuclear PER1 mRNA")
 MnPo annotate("units", "nM")

 val McPo = Double("c:McPo") default(0.00264)
 McPo annotate("description", "Cytoplasmic PER1 mRNA")
 McPo annotate("units", "nM")

 val trPo = Double("c:trPo") default(807.4)
 trPo annotate("description", "Transcription of PER1")
 trPo annotate("units", "1/h")

 val umPo = Double("c:umPo") default(6.21)
 umPo annotate("description", "Degradation of PER1 mRNA")
 umPo annotate("units", "1/h")


 /* PER2 Gene */

 val MnPt = Double("c:MnPt") default(0.015)
 MnPt annotate("description", "Nuclear PER2 mRNA")
 MnPt annotate("units", "nM")

 val McPt = Double("c:McPt") default(0.017)
 McPt annotate("description", "Cytoplasmic PER2 mRNA")
 McPt annotate("units", "nM")

 val trPt = Double("c:trPt") default(308.8)
 trPt annotate("description", "Transcription of PER2")
 trPt annotate("units", "1/h")

 val umPt = Double("c:umPt") default(0.38)
 umPt annotate("description", "Degradation of PER2 mRNA")
 umPt annotate("units", "1/h")


 /* PER Proteins */

 val Po = Double("c:Po") default(0.096)
 Po annotate("description", "Cytoplasmic PER1")
 Po annotate("units", "nM")

 val Pt = Double("c:Pt") default(1.018)
 Pt annotate("description", "Cytoplasmic PER2")
 Pt annotate("units", "nM")

 val PoC = Double("c:PoC") default(0.055)
 PoC annotate("description", "Cytoplasmic kinase-bound PER1")
 PoC annotate("units", "nM")

 val PtC = Double("c:PtC") default(0.604)
 PtC annotate("description", "Cytoplasmic kinase-bound PER2")
 PtC annotate("units", "nM")

 val tlp = Double("c:tlp") default(10)
 tlp annotate("description", "Translation of PER1 and PER2")
 tlp annotate("units", "1/h")

 val upu = Double("c:upu") default(0.08)
 upu annotate("description", "Degradation of unphosphorylated PER")
 upu annotate("units", "1/h")

 val PopC = Double("c:PopC") default(0.0012)
 PopC annotate("description", "Cytoplasmic kinase-bound PER1p")
 PopC annotate("units", "nM")

 val PtpC = Double("c:PtpC") default(0.005)
 PtpC annotate("description", "Cytoplasmic kinase-bound PER2p")
 PtpC annotate("units", "nM")

 val Pop = Double("c:Pop") default(0.000526)
 Pop annotate("description", "Cytoplasmic PER1p")
 Pop annotate("units", "nM")

 val Ptp = Double("c:Ptp") default(0.0027)
 Ptp annotate("description", "Cytoplasmic PER2p")
 Ptp annotate("units", "nM")

 val PonpCn = Double("c:PonpCn") default(0.0006)
 PonpCn annotate("description", "Nuclear kinase-bound PER1p")
 PonpCn annotate("units", "nM")

 val PtnpCn = Double("c:PtnpCn") default(0.0065)
 PtnpCn annotate("description", "Nuclear kinase-bound PER2p")
 PtnpCn annotate("units", "nM")

 val Ponp = Double("c:Ponp") default(0.0002)
 Ponp annotate("description", "Nuclear PER1p")
 Ponp annotate("units", "nM")

 val Ptnp = Double("c:Ptnp") default(0.0021)
 Ptnp annotate("descripion", "Nuclear PER2p")
 Ptnp annotate("units", "nM")

 val PoppC = Double("c:PoppC") default(0.0008)
 PoppC annotate("description", "Cytoplasmic kinase-bound PER1pp")
 PoppC annotate("units", "nM")

 val PtppC = Double("c:PtppC") default(0)
 PtppC annotate("description", "Cytoplasmic kinase-bound PER2pp")
 PtppC annotate("units", "nM")

 val Popp = Double("c:Popp") default(0.0011)
 Popp annotate("description", "Cytoplasmic PER1pp")
 Popp annotate("units", "nM")

 val Ptpp = Double("c:Ptpp") default(0)
 Ptpp annotate("description", "Cytoplasmic PER2pp")
 Ptpp annotate("units", "nM")

 val PonppCn = Double("c:PonppCn") default(0.000323)
 PonppCn annotate("description", "Nuclear kinase-bound PER1pp")
 PonppCn annotate("units", "nM")

 val PtnppCn = Double("c:PtnppCn") default(0)
 PtnppCn annotate("description", "Nuclear kinase-bound PER2pp")
 PtnppCn annotate("units", "nM")

 val Ponpp = Double("c:Ponpp") default(0.000088)
 Ponpp annotate("description", "Nuclear PER1pp")
 Ponpp annotate("units", "nM")

 val Ptnpp = Double("c:Ptnpp") default(0)
 Ptnpp annotate("description", "Nuclear PER2pp")
 Ptnpp annotate("units", "nM")


 /* CRY1 Gene */

 val MnRo = Double("c:MnRo") default(2.478)
 MnRo annotate("description", "Nuclear CRY1 mRNA")
 MnRo annotate("units", "nM")

 val McRo = Double("c:McRo") default(3.486)
 McRo annotate("description", "Cytoplasmic CRY1 mRNA")
 McRo annotate("units", "nM")

 val trRo = Double("c:trRo") default(9.03)
 trRt annotate("description", "Transcription of CRY1")
 trRt annotate("units", "1/h")


 /* CRY2 Gene */

 val MnRt = Double("c:MnRt") default(2.1)
 MnRt annotate("description", "Nuclear CRY2 mRNA")
 MnRt annotate("units", "nM")

 val McRt = Double("c:McRt") default(2.96)
 McRt annotate("description", "Cytoplasmic CRY2 mRNA")
 McRt annotate("units", "nM")

 val trRt = Double("c:trRt") default(7.66)
 trRt annotate("description", "Transcription of CRY2")
 trRt annotate("units", "1/h")


 /* CRY Proteins */

 val Ro = Double("c:Ro") default(8.2)
 Ro annotate("description", "Cytoplasmic CRY1")
 Ro annotate("units", "nM")

 val Rt = Double("c:Rt") default(5.18)
 Rt annotate("description", "Cytoplasmic CRY2")
 Rt annotate("units", "nM")

 val Ron = Double("c:Ron") default(0.0253)
 Ron annotate("description", "Nuclear CRY1")
 Ron annotate("units", "nM")

 val Rtn = Double("c:Rtn") default(0.0124)
 Rtn annotate("description", "Nuclear CRY2")
 Rtn annotate("units", "nM")

 val tlr = Double("c:tlr") default(1.031)
 tlr annotate("description", "Translation of CRY1 and CRY2")
 tlr annotate("units", "1/h")

 val uro = Double("c:uro") default(0.44)
 uro annotate("description", "Degradation of CRY1 unbound to PER")
 uro annotate("units", "1/h")

 val urt = Double("c:urt") default(0.59)
 urt annotate("description", "Degradation of CRY2 unbound to PER")
 urt annotate("units", "1/h")


 /* REVERBa Gene */

 val MnRv = Double("c:MnRv") default(0.000182)
 MnRv annotate("description", "Nuclear REVERBa mRNA")
 MnRv annotate("units", "nM")

 val McRv = Double("c:McRv") default(0.000005)
 McRv annotate("description", "Cytoplasmic REVERBa mRNA")
 McRv annotate("units", "nM")

 val trRv = Double("c:trRv") default(0.05)
 trRv annotate("description", "Transcription of REVERBa")
 trRv annotate("units", "1/h")

 val umRv = Double("c:umRv") default(15.11)
 umRv annotate("description", "Degradation of REVERBa mRNA")
 umRv annotate("units", "1/h")


 /* REVERBa Protein */

 val Rv = Double("c:Rv") default(0.000001)
 Rv annotate("description", "Cytoplasmic REVERBa")
 Rv annotate("units", "nM")

 val Rvn = Double("c:Rvn") default(0)
 Rvn annotate("description", "Nuclear REVERBa")
 Rvn annotate("units", "nM")

 val RvRv = Double("c:RvRv") default(0)
 RvRv annotate("description", "Cytoplasmic REVERBa dimer")
 RvRv annotate("units", "nM")

 val RvnRvn = Double("c:RvnRvn") default(0)
 RvnRvn annotate("description", "Nuclear REVERBa dimer")
 RvnRvn annotate("units", "nM")

 val tlrv = Double("c:tlrv") default(2.53)
 tlrv annotate("description", "Translation of REVERBa")
 tlrv annotate("units", "1/h")

 val arv = Double("c:arv") default(0.21)
 arv annotate("description", "Dimerisation of REVERBa")
 arv annotate("units", "1/nM 1/h")

 val drv = Double("c:drv") default(3.62)
 drv annotate("description", "Undimerisation of REVERBa")
 drv annotate("units", "1/h")

 val uRv = Double("c:uRv") default(16.25)
 uRv annotate("description", "Degradation of REVERBa")


 /* PER-CRY Complexes */

 val PopRo = Double("c:PopRo") default(0.0009)
 PopRo annotate("description", "Cytoplasmic PER1p-CRY1")
 PopRo annotate("units", "nM")

 val PtpRo = Double("c:PtpRo") default(0.01)
 PtpRo annotate("description", "Cytoplasmic PER2p-CRY1")
 PtpRo annotate("units", "nM")

 val PopCRo = Double("c:PopCRo") default(0.0007)
 PopCRo annotate("description", "Cytoplasmic kinase-bound PER1p-CRY1")
 PopCRo annotate("units", "nM")

 val PtpCRo = Double("c:PtpCRo") default(0.009)
 PtpCRo annotate("description", "Cytoplasmic kinase-bound PER2p-CRY1")
 PtpCRo annotate("units", "nM")

 val PonpRon = Double("c:PonpRon") default(0.000543)
 PonpRon annotate("description", "Nuclear PER1p-CRY1")
 PonpRon annotate("units", "nM")

 val PtnpRon = Double("c:PtnpRon") default(0.012)
 PtnpRon annotate("description", "Nuclear PER2p-CRY1")
 PtnpRon annotate("units", "nM")

 val PonpCnRon = Double("c:PonpCnRon") default(0.0015)
 PonpCnRon annotate("description", "Nuclear kinase-bound PER1p-CRY1")
 PonpCnRon annotate("units", "nM")

 val PtnpCnRon = Double("c:PtnpCnRon") default(0.0373)
 PtnpCnRon annotate("description", "Nuclear kinase-bound PER2p-CRY1")
 PtnpCnRon annotate("units", "nM")

 val PopRt = Double("c:PopRt") default(0.00055)
 PopRt annotate("description", "Cytoplasmic PER1p-CRY2")
 PopRt annotate("units", "nM")

 val PtpRt = Double("c:PtpRt") default(0.0064)
 PtpRt annotate("description", "Cytoplasmic PER2p-CRY2")
 PtpRt annotate("units", "nM")

 val PopCRt = Double("c:PopCRt") default(0.0004)
 PopCRt annotate("description", "Cytoplasmic kinase-bound PER1p-CRY2")
 PopCRt annotate("units", "nM")

 val PtpCRt = Double("c:PtpCRt") default(0.0055)
 PtpCRt annotate("description", "Cytoplasmic kinase-bound PER2p-CRY2")
 PtpCRt annotate("units", "nM")

 val PonpRtn = Double("c:PonpRtn") default(0.00033)
 PonpRtn annotate("description", "Nuclear PER1p-CRY2")
 PonpRtn annotate("units", "nM")

 val PtnpRtn = Double("c:PtnpRtn") default(0.007)
 PtnpRtn annotate("description", "Nuclear PER2p-CRY2")
 PtnpRtn annotate("units", "nM")

 val PonpCnRtn = Double("c:PonpCnRtn") default(0.0009)
 PonpCnRtn annotate("description", "Nuclear kinase-bound PER1p-CRY2")
 PonpCnRtn annotate("units", "nM")

 val PtnpCnRtn = Double("c:PtnpCnRtn") default(0.022)
 PtnpCnRtn annotate("description", "Nuclear kinase-bound PER2p-CRY2")
 PtnpCnRtn annotate("units", "nM")

 val PoppRo = Double("c:PoppRo") default(0.0124)
 PoppRo annotate("description", "Cytoplasmic PER1pp-CRY1")
 PoppRo annotate("units", "nM")

 val PtppRo = Double("c:PtppRo") default(0)
 PtppRo annotate("description", "Cytoplasmic PER2pp-CRY1")
 PtppRo annotate("units", "nM")

 val PoppCRo = Double("c:PoppCRo") default(0.008)
 PoppCRo annotate("description", "Cytoplasmic kinase-bound PER1pp-CRY1")
 PoppCRo annotate("units", "nM")

 val PtppCRo = Double("c:PtppCRo") default(0)
 PtppCRo annotate("description", "Cytoplasmic kinase-bound PER2pp-CRY1")
 PtppCRo annotate("units", "nM")

 val PonppRon = Double("c:PonppRon") default(0.000525)
 PonppRon annotate("description", "Nuclear PER1pp-CRY1")
 PonppRon annotate("units", "nM")

 val PtnppRon = Double("c:PtnppRon") default(0)
 PtnppRon annotate("description", "Nuclear PER2pp-CRY1")
 PtnppRon annotate("units", "nM")

 val PonppCnRon = Double("c:PonppCnRon") default(0.0018)
 PonppCnRon annotate("description", "Nuclear kinase-bound PER1pp-CRY1")
 PonppCnRon annotate("units", "nM")

 val PtnppCnRon = Double("c:PtnppCnRon") default(0)
 PtnppCnRon annotate("description", "Nuclear kinase-bound PER2pp-CRY1")
 PtnppCnRon annotate("units", "nM")

 val PoppRt = Double("c:PoppRt") default(0.0077)
 PoppRt annotate("description", "Cytoplasmic PER1pp-CRY2")
 PoppRt annotate("units", "nM")

 val PtppRt = Double("c:PtppRt") default(0)
 PtppRt annotate("description", "Cytoplasmic PER2pp-CRY2")
 PtppRt annotate("units", "nM")

 val PoppCRt = Double("c:PoppCRt") default(0.005)
 PoppCRt annotate("description", "Cytoplasmic PER1pp-CRY2")
 PoppCRt annotate("units", "nM")

 val PtppCRt = Double("c:PtppCRt") default(0)
 PtppCRt annotate("description", "Cytoplasmic PER2pp-CRY2")
 PtppCRt annotate("units", "nM")

 val PonppRtn = Double("c:PonppRtn") default(0.00031)
 PonppRtn annotate("description", "Nuclear PER1pp-CRY2")
 PonppRtn annotate("units", "nM")

 val PtnppRtn = Double("c:PtnppRtn") default(0)
 PtnppRtn annotate("description", "Nuclear PER2pp-CRY2")
 PtnppRtn annotate("units", "nM")

 val PonppCnRtn = Double("c:PonppCnRtn") default(0.001)
 PonppCnRtn annotate("description", "Nuclear kinase-bound PER1pp-CRY2")
 PonppCnRtn annotate("units", "nM")

 val PtnppCnRtn = Double("c:PtnppCnRtn") default(0)
 PtnppCnRtn annotate("description", "Nuclear kinase-bound PER2pp-CRY2")
 PtnppCnRtn annotate("units", "nM")



 d(G) := (bin * Rn * (1-G)) - (unbin * G)
 d(GRv) := (binRv * RvnRvn * (1-GRv)) - (unbinRv * GRv)

 d(MnRo) := (trRo * (1-G) * (1-GRv)^3) - (tmc * MnRo)
 d(McRo) := (tmc * MnRo) - (umR * McRo)
 d(MnRt) := (trRt * (1-G)) - (tmc * MnRt)
 d(McRt) := (tmc * MnRt) - (umR * McRt)
 d(MnPo) := (trPo * (1-G)^5) + L) - (tmc * MnPo)
 d(McPo) := (tmc * MnPo) - (umPo * McPo)
 d(MnPt) := (trPt * (1-G)^5) + L) - (tmc * MnPt)
 d(McPt) := (tmc * MnPt) - (umPt * McPt)
 d(MnRv) := (trRv * (1-G)^3) - (tmc * MnRv)
 d(McRv) := (tmc * MnRv) - (umRv * McRv)

 d(Rv) := (tlrv * McRv) - (2 * arv * Rv * Rv) + (2 * drv * RvRv) - (nl * Rv) + (ne * Rvn) - (uRv * Rv)
 d(Rvn) := - (2 * Nf * arv * Rvn * Rvn) + (2 * drv * RvnRvn) + (nl * Rv) - (ne * Rvn) - (uRv * Rvn)
 d(RvRv) := (arv * Rv * Rv) - (drv * RvRv) - (nl * RvRv) + (ne * RvnRvn) - (2 * uRv * RvRv)
 d(RvnRvn) := (Nf * arv * Rvn * Rvn) - (drv * RvnRvn) + (nl * RvRv) - (ne * RvnRvn) -
(2 * uRv * RvnRvn)

 d(Po) := (tlp * McPo) - (ac * Po * C) + (dc * PoC) - (upu * Po)
 d(Pt) := (tlp * McPt) - (ac * Pt * C) + (dc * PtC) - (upu * Pt)
 d(PoC) := (ac * Po * C) - (dc * PoC) - (hoo * PoC) - (upu * PoC)
 d(PtC) := (ac * Pt * C) - (dc * PtC) - (hot * PtC) - (upu * PtC)
 d(PopC) := (hoo * PoC) + (ac * Pop * C) - (dc * PopC) - (up * PopC) - (hto * PopC) - (nl * PopC) + (ne * PonpCn) - (ar * PopC * Ro) + (dr * PopCRo) - (ar * PopC * Rt) + (dr * PopCRt)
 d(PtpC) := (hot * PtC) + (ac * Ptp * C) - (dc * PtpC) - (up * PtpC) - (ht * PtpC) - (nl * PtpC) + (ne * PtnpCn) - (ar * PtpC * Ro) + (dr * PtpCRo) - (ar * PtpC * Rt) + (dr * PtpCRt)
 d(Pop) := - (ac * Pop * C) + (dc * PopC) - (up * Pop) - (ar * Pop * Ro) + (dr * PopRo) - (ar * Pop * Rt) + (dr * PopRt) - (nl * Pop) + (ne * Ponp)
 d(Ptp) := - (ac * Ptp * C) + (dc * PtpC) - (up * Ptp) - (ar * Ptp * Ro) + (dr * PtpRo) - (ar * Ptp * Rt) + (dr * PtpRt) - (nl * Ptp) + (ne * Ptnp)
 d(PoppC) := (hto * PopC) - (up * PoppC) + (ac * Popp * C) - (dc * PoppC) + (ne * PonppCn) - (ar * PoppC * Ro) + (dr * PoppCRo) - (ar * PoppC * Rt) + (dr * PoppCRt)
 d(PtppC) := (ht * PtpC) - (up * PtppC) + (ac * Ptpp * C) - (dc * PtppC) + (ne * PtnppCn) - (ar * PtppC * Ro) + (dr * PtppCRo) - (ar * PtppC * Rt) + (dr * PtppCRt)
 d(PopRo) := (ar * Pop * Ro) - (dr * PopRo) - (ac * PopRo * C) + (dc * PopCRo) - (nl * PopRo) + (ne * PonpRon)
 d(PtpRo) := (ar * Ptp * Ro) - (dr * PtpRo) - (ac * PtpRo * C) + (dc * PtpCRo) - (nl * PtpRo) + (ne * PtnpRon)
 d(PopRt) := (ar * Pop * Rt) - (dr * PopRt) - (ac * PopRt * C) + (dc * PopCRt) - (nl * PopRt) + (ne * PonpRtn)
 d(PtpRt) := (ar * Ptp * Rt) - (dr * PtpRt) - (ac * PtpRt * C) + (dc * PtpCRt) - (nl * PtpRt) + (ne * PtnpRtn)
 d(PoppRo) := (ar * Popp * Ro) - (dr * PoppRo) - (ac * PoppRo * C) + (dc * PoppCRo) + (ne * PonppRon)
 d(PoppRt) := (ar * Popp * Rt) - (dr * PoppRt) - (ac * PoppRt * C) + (dc * PoppCRt) + (ne * PonppRtn)
 d(PtppRo) := (ar * Ptpp * Ro) - (dr * PtppRo) - (ac * PtppRo * C) + (dc * PtppCRo) + (ne * PtnppRon)
 d(PtppRt) := (ar * Ptpp * Rt) - (dr * PtppRt) - (ac * PtppRt * C) + (dc * PtppCRt) + (ne * PtnppRtn)
 d(Popp) := - (ac * Popp * C) + (dc * PoppC) + (ne * Ponpp) - (ar * Popp * Ro) + (dr * PoppRo) – (ar * Popp * Rt) + (dr * PoppRt) - (up * Popp)
 d(Ptpp) := - (ac * Ptpp * C) + (dc * PtppC) + (ne * Ptnpp) - (ar * Ptpp * Ro) + (dr * PtppRo) - (ar * Ptpp * Rt) + (dr * PtppRt) - (up * Ptpp)
 d(PopCRo) := (ar * PopC * Ro) - (dr * PopCRo) + (ac * PopRo * C) - (dc * PopCRo) - (nl * PopCRo) + (ne * PonpCnRon) - (hto * PopCRo)
 d(PtpCRo) := (ar * PtpC * Ro) - (dr * PtpCRo) + (ac * PtpRo * C) - (dc * PtpCRo) - (nl * PtpCRo) + (ne * PtnpCnRon) - (ht * PtpCRo)
 d(PopCRt) := (ar * PopC * Rt) - (dr * PopCRt) + (ac * PopRt * C) - (dc * PopCRt) - (nl * PopCRt) + (ne * PonpCnRtn) - (hto * PopCRt)
 d(PtpCRt) := (ar * PtpC * Rt) - (dr * PtpCRt) + (ac * PtpRt * C) - (dc * PtpCRt) - (nl * PtpCRt) + (ne * PtnpCnRtn) - (ht * PtpCRt)
 d(PoppCRo) := (ar * PoppC * Ro) - (dr * PoppCRo) + (ac * PoppRo * C) - (dc * PoppCRo) + (ne * PonppCnRon) + (hto * PopCRo)
 d(PtppCRo) := (ar * PtppC * Ro) - (dr * PtppCRo) + (ac * PtppRo * C) - (dc * PtppCRo) + (ne * PtnppCnRon) + (ht * PtpCRo)
 d(PoppCRt) := (ar * PoppC * Rt) - (dr * PoppCRt) + (ac * PoppRt * C) - (dc * PoppCRt) + (ne * PonppCnRtn) + (hot * PopCRt)
 d(PtppCRt) := (ar * PtppC * Rt) - (dr * PtppCRt) + (ac * PtppRt * C) - (dc * PtppCRt) + (ne * PtnppCnRtn) + (ht * PtpCRt)

 d(Ro) := - (ar * Ro * Pop) - (ar * Ro * Popp) - (ar * Ro * PopC) - (ar * Ro * PoppC) + (dr * PopRo) + (dr * PoppRo) + (dr * PopCRo) + (dr * PoppCRo) - (ar * Ro * Ptp) - (ar * Ro * Ptpp) - (ar * Ro * PtpC) - (ar * Ro * PtppC) + (dr * PtpRo) + (dr * PtppRo) + (dr * PtpCRo) + (dr * PtppCRo) + (tlr * McRo) - (uro * Ro)
 d(Rt) := - (ar * Rt * Pop) - (ar * Rt * Popp) - (ar * Rt * PopC) - (ar * Rt * PoppC) + (dr * PopRt) + (dr * PoppRt) + (dr * PopCRt) + (dr * PoppCRt) - (ar * Rt * Ptp) - (ar * Rt * Ptpp) - (ar * Rt * PtpC) - (ar * Rt * PtppC) + (dr * PtpRt) + (dr * PtppRt) + (dr * PtpCRt) + (dr * PtppCRt) + (tlr * McRt) - (urt * Rt)

 d(PonpCn) := (ac * Nf * Ponp * Cn) - (dc * PonpCn) - (hto * PonpCn) + (nl * PopC) - (ne * PonpCn) - (ar * Nf * PonpCn * Ron) + (dr * PonpCnRon) - (ar * Nf * PonpCn * Rtn) + (dr * PonpCnRtn) - (up * PonpCn)
 d(PtnpCn) := (ac * Nf * Ptnp * Cn) - (dc * PtnpCn) - (ht * PtnpCn) + (nl * PtpC) - (ne * PtnpCn) - (ar * Nf * PtnpCn * Ron) + (dr * PtnpCnRon) - (ar * Nf * PtnpCn * Rtn) + (dr * PtnpCnRtn) - (up * PtnpCn)
 d(Ponp) := - (ac * Nf * Ponp * Cn) + (dc * PonpCn) - (ar * Nf * Ponp * Ron) + (dr * PonpRon) - (ar * Nf * Ponp * Rtn) + (dr * PonpRtn) + (nl * Pop) - (ne * Ponp) - (up * Ponp)
 d(Ptnp) := - (ac * Nf * Ptnp * Cn) + (dc * PtnpCn) - (ar * Nf * Ptnp * Ron) + (dr * PtnpRon) - (ar * Nf * Ptnp * Rtn) + (dr * PtnpRtn) + (nl * Ptp) - (ne * Ptnp) - (up * Ptnp)
 d(PonppCn) := (hto * PonpCn) + (ac * Nf * Ponpp * Cn) - (dc * PonppCn) - (ne * PonppCn) - (ar * Nf * PonppCn * Ron) + (dr * PonppCnRon) – (ar * Nf * PonppCn * Rtn) + (dr * PonppCnRtn) - (up * PonppCn)
 d(PtnppCn) := (ht * PtnpCn) + (ac * Nf * Ptnpp * Cn) - (dc * PtnppCn) - (ne * PtnppCn) - (ar * Nf * PtnppCn * Ron) + (dr * PtnppCnRon) - (ar * Nf * PtnppCn * Rtn) + (dr * PtnppCnRtn) - (up * PtnppCn)
 d(PonpRon) := (ar * Nf * Ponp * Ron) - (dr * PonpRon) - (ac * Nf * PonpRon * Cn) + (dc * PonpCnRon) + (nl * PopRo) - (ne * PonpRon)
 d(PtnpRon) := (ar * Nf * Ptnp * Ron) - (dr * PtnpRon) - (ac * Nf * PtnpRon * Cn) + (dc * PtnpCnRon) + (nl * PtpRo) - (ne * PtnpRon)
 d(PonpRtn) := (ar * Nf * Ponp * Rtn) - (dr * PonpRtn) - (ac * Nf * PonpRtn * Cn) + (dc * PonpCnRtn) + (nl * PopRt) - (ne * PonpRtn)
 d(PtnpRtn) := (ar * Nf * Ptnp * Rtn) - (dr * PtnpRtn) - (ac * Nf * PtnpRtn * Cn) + (dc * PtnpCnRtn) + (nl * PtpRt) - (ne * PtnpRtn)
 d(PonppRon) := (ar * Nf * Ponpp * Ron) - (dr * PonppRon) - (ac * Nf * PonppRon * Cn) + (dc * PonppCnRon) - (ne * PonppRon)
 d(PtnppRon) := (ar * Nf * Ptnpp * Ron) - (dr * PtnppRon) - (ac * Nf * PtnppRon * Cn) + (dc * PtnppCnRon) - (ne * PtnppRon)
 d(PonppRtn) := (ar * Nf * Ponpp * Rtn) - (dr * PonppRtn) - (ac * Nf * PonppRtn * Cn) + (dc * PonppCnRtn) - (ne * PonppRtn)
 d(PtnppRtn) := (ar * Nf * Ptnpp * Rtn) - (dr * PtnppRtn) - (ac * Nf * PtnppRtn * Cn) + (dc * PtnppCnRtn) - (ne * PtnppRtn)
 d(Ponpp) := - (ac * Nf * Ponpp * Cn) + (dc * PonppCn) - (ne * Ponpp) - (ar * Nf * Ponpp * Ron) + (dr * PonppRon) - (ar * Nf * Ponpp * Rtn) + (dr * PonppRtn) - (up * Ponpp)
 d(Ptnpp) := - (ac * Nf * Ptnpp * Cn) + (dc * PtnppCn) - (ne * Ptnpp) - (ar * Nf * Ptnpp * Ron) + (dr * PtnppRon) - (ar * Nf * Ptnpp * Rtn) + (dr * PtnppRtn) - (up * Ptnpp)
 d(PonpCnRon) := (ar * Nf * PonpCn * Ron) - (dr * PonpCnRon) + (ac * Nf * PonpRon * Cn) - (dc * PonpCnRon) + (nl * PopCRo) - (ne * PonpCnRon) - (hto * PonpCnRon)
 d(PtnpCnRon) := (ar * Nf * PtnpCn * Ron) - (dr * PtnpCnRon) + (ac * Nf * PtnpRon * Cn) - (dc * PtnpCnRon) + (nl * PtpCRo) - (ne * PtnpCnRon) - (ht * PtnpCnRon)
 d(PonpCnRtn) := (ar * Nf * PonpCn * Rtn) - (dr * PonpCnRtn) + (ac * Nf * PonpRtn * Cn) - (dc * PonpCnRtn) + (nl * PopCRt) - (ne * PonpCnRtn) - (hto * PonpCnRtn)
 d(PtnpCnRtn) := (ar * Nf * PtnpCn * Rtn) - (dr * PtnpCnRtn) + (ac * Nf * PtnpRtn * Cn) - (dc * PtnpCnRtn) + (nl * PtpCRt) - (ne * PtnpCnRtn) - (ht * PtnpCnRtn)
 d(PonppCnRon) := (ar * Nf * PonppCn * Ron) - (dr * PonppCnRon) + (ac * Nf * PonppRon * Cn) - (dc * PonppCnRon) - (ne * PonppCnRon) + (hto * PonpCnRon)
 d(PtnppCnRon) := (ar * Nf * PtnppCn * Ron) - (dr * PtnppCnRon) + (ac * Nf * PtnppRon * Cn) - (dc * PtnppCnRon) - (ne * PtnppCnRon) + (ht * PtnpCnRon)
 d(PonppCnRtn) := (ar * Nf * PonppCn * Rtn) - (dr * PonppCnRtn) + (ac * Nf * PonppRtn * Cn) - (dc * PonppCnRtn) - (ne * PonppCnRtn) + (hto * PonpCnRtn)
 d(PtnppCnRtn) := (ar * Nf * PtnppCn * Rtn) - (dr * PtnppCnRtn) + (ac * Nf * PtnppRtn * Cn) - (dc * PtnppCnRtn) - (ne * PtnppCnRtn) + (ht * PtnpCnRtn)

 d(Ron) := - (ar * Nf * Ron * Ponp) - (ar * Nf * Ron * Ponpp) - (ar * Nf * Ron * PonpCn) - (ar * Nf * Ron * PonppCn) + (dr * PonpRon) + (dr * PonppRon) + (dr * PonpCnRon) + (dr * PonppCnRon) - (ar * Nf * Ron * Ptnp) - (ar * Nf * Ron * Ptnpp) - (ar * Nf * Ron * PtnpCn) - (ar * Nf * Ron * PtnppCn) + (dr * PtnpRon) + (dr * PtnppRon) + (dr * PtnpCnRon) + (dr * PtnppCnRon) - (uro * Ron)
 d(Rtn) := - (ar * Nf * Rtn * Ponp) - (ar * Nf * Rtn * Ponpp) - (ar * Nf * Rtn * PonpCn) - (ar * Nf * Rtn * PonppCn) + (dr * PonpRtn) + (dr * PonppRtn) + (dr * PonpCnRtn) + (dr * PonppCnRtn) - (ar * Nf * Rtn * Ptnp) - (dr * Nf * Rtn * Ptnpp) - (ar * Nf * Rtn * PtnpCn) - (ar * Nf * Rtn * PtnppCn) + (dr * PtnpRtn) + (dr * PtnppRtn) + (dr * PtnpCnRtn) + (dr * PtnppCnRtn) - (urt * Rtn)

 d(Cn) := - (ac * Nf * Cn * Ponp) - (ac * Nf * Cn * Ponpp) - (ac * Nf * Cn * PonpRon) - (ac * Nf * Cn * PonppRon) + (dc * PonpCn) + (dc * PonppCn) + (dc * PonpCnRon) + (dc * PonppCnRon) - (ac * Nf * Cn * Ptnp) - (ac * Nf * Cn * Ptnpp) - (ac * Nf * Cn * PtnpRon) - (ac * Nf * Cn * PtnppRon) + (dc * PtnpCn) + (dc * PtnppCn) + (dc * PtnpCnRon) + (dc * PtnppCnRon) - (ac * Nf * Cn * PonpRtn) - (ac * Nf * Cn * PonppRtn) + (dc * PonpCnRtn) + (dc * PonppCnRtn) - (ac * Nf * Cn * PtnpRtn) - (ac * Nf * Cn * PtnppRtn) + (dc * PtnpCnRtn) + (dc * PtnppCnRtn) + (up * PonpCn) + (up * PonppCn) + (up * PtnpCn) + (up * PtnppCn)

}


class MammalianCircadianClockModel extends Model {

 /* Handle light level at this stage, feeding it to process?
  * Can we access simulation time for perturbations?
  * 
  *  val Lon = Double("c:Lon") default(0.000339)
  *  Lon annotate("description", "Light level (light on)")
  *
  *  val Loff = Double("c:Loff") default(0)
  *  Loff annotate("description", "Light level (light off)")
  *
  * If Math.floor([T]/12) is odd, use Loff as input to MammalianCircadianClock.
  * Else, use Lon.
  *
  */

 val process = new MammalianCircadianClock(L)

  /* Both quantities below should remain constant during simulation:
   *
   *  C = Ct - (PoC + PtC + PopC + PtpC + PoppC + PtppC + PopCRo + PopCRt + PtpCRo + PtpCRt + PoppCRo + PoppCRt + PtppCRo + PtppCRt + PonpCn + PtnpCn + PonppCn
+ PtnppCn + PonpCnRon + PonpCnRtn + PtnpCnRon + PtnpCnRtn + PonppCnRon + PonppCnRtn + PtnppCnRon + PtnppCnRtn + Cn)
   *
   *  Rn = (Ron + PonpRon + PonppRon + PonpCnRon + PonppCnRon + PtnpRon + PtnppRon + PtnpCnRon + PtnppCnRon + Rtn + PonpRtn + PonppRtn + PonpCnRtn + PonppCnRtn + PtnpRtn + PtnppRtn + PtnpCnRtn + PtnppCnRtn
   *
   */

}
