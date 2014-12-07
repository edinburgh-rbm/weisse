/*
 *  CellModel by Weisse et al.
 *  Copyright (C) 2014  Andrea Y. Weisse
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package uk.ac.ed.inf.mois.weisse

import uk.ac.ed.inf.mois.{Model, Process, ProcessGroup}
import uk.ac.ed.inf.mois.sched.CompositionScheduler
import uk.ac.ed.inf.mois.{VarCalc, Math}
import uk.ac.ed.inf.mois.reaction.DeterministicReactionNetwork
import spire.implicits._
import uk.ac.ed.inf.mois.implicits._


class JWeisseRates(val k_cm : Double, val Kp : Double, val Kt : Double, val Km : Double, val M : Double, val gmax : Double, val nr : Double, val nx : Double, val cl : Double, val s0 : Double, val vm : Double, val vt : Double
) extends VarCalc {
  /* ATP and internal nutrient */
  val a = Double("a") default(1.0)
  val si = Double("si")

  /* proteins */
  val r = Double("r") default(1.0)
  val p = Double("p")
  val q = Double("q")
  val em = Double("em")
  val et = Double("et")

  /* ribosome-bound mRNA */
  val rmr = Double("rmr")
  val rmt = Double("rmt")
  val rmm = Double("rmm")
  val rmp = Double("rmp")
  val rmq = Double("rmq")

  /* ribosome-bound mRNA sequestered by chloramphenicol */
  val zmr = Double("zmr")
  val zmt = Double("zmt")
  val zmm = Double("zmm")
  val zmp = Double("zmp")
  val zmq = Double("zmq")

  /* some rates */
  val Kgamma = Double("Kgamma") nonnegative()
  val gamma = Double("gamma") nonnegative()
  val ttrate = Double("ttrate") nonnegative()
  val lam = Double("lam") nonnegative()
  val fr = Double("fr") nonnegative()
  val nucat = Double("nucat") nonnegative()
  val nurat = Double("nurat") nonnegative()
  val f = Double("f") nonnegative()
  val b = Double("b") nonnegative()

  calc(Kgamma) := gmax/Kp
  calc(gamma) := gmax*a/(Kgamma + a)
  calc(ttrate) := gamma*(rmq + rmr + rmp + rmt + rmm)
  calc(lam) := ttrate/M
  calc(fr) := nr*(r + rmr + rmp + rmt + rmm + rmq + zmr + zmp + zmt + zmm + zmq) /
    ( nr*(r + rmr + rmp + rmt + rmm + rmq + zmr + zmp + zmt + zmm + zmq) + nx*(p + q + et + em) )
  calc(nucat) := em*vm*si/(Km + si)
  calc(nurat) := et*vt*s0/(Kt + s0)
  calc(f) := cl*k_cm
  calc(b) := 0.0
}


class JWeisseModel extends Model {

  /* parameters from MCMC sampling */
  val Kp = 180.1378030928276
  val thetar = 426.8693338968694
  val k_cm = 0.005990373118888
  val wr = 929.9678874564831
  val wq = 948.9349882947897
  val we = 4.139172187824451
  val Kq = 1.522190403737490e+05
  val thetax = 4.379733394834643

  /* no p-transcription (from Copasi's deterministic optimization) */
  val wp = 0.0

  /* other parameters mostly from literature */
  val gmax = 1260.0
  val M = 1.0e8
  val vt = 726.0
  val Kt = 1.0e3
  val s0 = 1.0e4
  val vm = 5800.0
  val Km = 1.0e3
  val ns = 0.5
  val nq = 4.0
  val nr = 7549.0
  val nx = 300.0
  val dm = 0.1
  val cl = 0.0
  val kb = 1.0
  val ku = 1.0

  val process = new ProcessGroup {
    scheduler = new CompositionScheduler(0.01)
  }

  process += new JWeisseRates(k_cm, Kp, Kt, Km, M, gmax, nr, nx, cl, s0, vm, vt)
  process += new JWeisseCellTranscription(thetar, thetax, wr, wq, we, wp, Kq, nq)
  process += new JWeisseCellTranslation(nr, nx)
  process += new JWeisseCellmRNADegradation(dm)
  process += new JWeisseCellmRNARibosomeComplex(kb, ku)
  process += new JWeisseCellChloramphenicol
  process += new JWeisseCellDilution
  process += new JWeisseCellMetabolism(ns)

  // -- Initial values -- (not communicating downwards - need to write defaults in process proper)

  Double("a") default(1.0)
  Double("r") default(1.0)

  // Everything else is just zero
}
