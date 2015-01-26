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
import uk.ac.ed.inf.mois.ode.{Rosenbrock, ODEDebugHandler}
import uk.ac.ed.inf.mois.sched.CompositionScheduler
import uk.ac.ed.inf.mois.{VarCalc, Math}
import uk.ac.ed.inf.mois.reaction.DeterministicReactionNetwork
import spire.math.Jet
import spire.implicits._
import uk.ac.ed.inf.mois.implicits._

class WeisseRatesSteady(
  val k_cm  : Double,
  val Kp    : Double,
  val Kt    : Double,
  val Km    : Double,
  val M     : Double,
  val gmax  : Double,
  val nr    : Double,
  val nx    : Double,
  val cl    : Double,
  val s0    : Double,
  val vm    : Double,
  val vt    : Double
) extends VarCalc {
  /* ATP and internal nutrient */
  val a = Double("a")
  val si = Double("si")

  /* proteins */
  val r = Double("r")
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

class WeisseCellSteady(
  val Kp    : Double,
  val thetar: Double,
  val k_cm  : Double,
  val wr    : Double,
  val wq    : Double,
  val we    : Double,
  val Kq    : Double,
  val thetax: Double,
  val wp    : Double,
  val gmax  : Double,
  val M     : Double,
  val vt    : Double,
  val Kt    : Double,
  val s0    : Double,
  val vm    : Double,
  val Km    : Double,
  val ns    : Double,
  val nq    : Double,
  val nr    : Double,
  val nx    : Double,
  val dm    : Double,
  val cl    : Double,
  val kb    : Double,
  val ku    : Double)
    extends DeterministicReactionNetwork[Double, Jet[Double]] with Rosenbrock
       with Math {

  /* define variables */
  /* ATP and internal nutrient */
  val a = Species("a") nonnegative()
  val si = Species("si") nonnegative()

  /* proteins */
  val r = Species("r") nonnegative()
  val et = Species("et") nonnegative()
  val em = Species("em") nonnegative()
  val p = Species("p") nonnegative()
  val q = Species("q") nonnegative()

  /* mRNA */
  val mr = Species("mr") nonnegative()
  val mt = Species("mt") nonnegative()
  val mm = Species("mm") nonnegative()
  val mp = Species("mp") nonnegative()
  val mq = Species("mq") nonnegative()

  /* ribosome-bound mRNA */
  val rmr = Species("rmr") nonnegative()
  val rmt = Species("rmt") nonnegative()
  val rmm = Species("rmm") nonnegative()
  val rmp = Species("rmp") nonnegative()
  val rmq = Species("rmq") nonnegative()

  /* ribosome-bound mRNA sequestered by chloramphenicol */
  val zmr = Species("zmr") nonnegative()
  val zmt = Species("zmt") nonnegative()
  val zmm = Species("zmm") nonnegative()
  val zmp = Species("zmp") nonnegative()
  val zmq = Species("zmq") nonnegative()

  /* some rates */
  val Kgamma = Double("Kgamma")
  val gamma = Double("gamma")
  val ttrate = Double("ttrate")
  val lam = Double("lam")
  val fr = Double("fr")
  val nucat = Double("nucat")
  val nurat = Double("nurat")
  val f = Double("f")
  val b = Double("b")

  import spire.math.Jet
  implicit def jetToDouble(j: Jet[Double]): Double = j.real

  reactions(
    /* nutrient import */
    () --> si `at!` nurat,

    /* nutrient metabolism */
    si --> () `at!` nucat,
    () --> a `at!` ns*nucat,
    a --> () `at!` ttrate,

    /* transcription */
    () --> mr `at!` wr*a/(thetar + a),
    () --> mt `at!` we*a/(thetax + a),
    () --> mm `at!` we*a/(thetax + a),
    () --> mp `at!` wp*a/(thetax + a),
    () --> mq `at!` wq*a/(thetax + a)/(1 + pow((q/Kq), nq)),

    /* translation */
    r + mr --> rmr at kb,
    r + mt --> rmt at kb,
    r + mm --> rmm at kb,
    r + mp --> rmp at kb,
    r + mq --> rmq at kb,

    rmr --> r + mr at ku,
    rmt --> r + mt at ku,
    rmm --> r + mm at ku,
    rmp --> r + mp at ku,
    rmq --> r + mq at ku,

    rmr --> r + r + mr `at!` gamma/nr,
    rmt --> r + et + mt `at!` gamma/nx,
    rmm --> r + em + mm `at!` gamma/nx,
    rmp --> r + p + mp `at!` gamma/nx,
    rmq --> r + q + mq `at!` gamma/nx,

    /* mRNA degradation */
    mr --> () at dm,
    mt --> () at dm,
    mm --> () at dm,
    mp --> () at dm,
    mq --> () at dm,

    /* chloramphenicol */
    rmr --> zmr at f,
    rmt --> zmt at f,
    rmm --> zmm at f,
    rmp --> zmp at f,
    rmq --> zmq at f,

    zmr --> rmr at b,
    zmt --> rmt at b,
    zmm --> rmm at b,
    zmp --> rmp at b,
    zmq --> rmq at b,

    /* dilution */
    mr --> () at lam,
    mt --> () at lam,
    mm --> () at lam,
    mp --> () at lam,
    mq --> () at lam,
    rmr --> () at lam,
    rmt --> () at lam,
    rmm --> () at lam,
    rmp --> () at lam,
    rmq --> () at lam,
    r --> () at lam,
    et --> () at lam,
    em --> () at lam,
    p --> () at lam,
    q --> () at lam,
    a --> () at lam,
    si --> () at lam,
    zmr --> () at lam,
    zmt --> () at lam,
    zmm --> () at lam,
    zmp --> () at lam,
    zmq --> () at lam
  )
}

class WeisseModelSteady extends Model {

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
    scheduler = new CompositionScheduler(1.0)
  }

  process += new WeisseCellSteady(
      Kp,
      thetar,
      k_cm,
      wr,
      wq,
      we,
      Kq,
      thetax,
      wp,
      gmax,
      M,
      vt,
      Kt,
      s0,
      vm,
      Km,
      ns,
      nq,
      nr,
      nx,
      dm,
      cl,
      kb,
      ku
  )

  process += new WeisseRatesSteady(
    k_cm, Kp, Kt, Km, M, gmax, nr, nx, cl, s0, vm, vt
  )

  // -- Initial values --

  /* ATP and internal nutrient */
  Double("a") default(43297.502)
  Double("si") default(31096.192)

  /* proteins */
  Double("r") default(0.0414)
  Double("et") default(471.364)
  Double("em") default(471.364)
  Double("p") default(471.364)
  Double("q") default(471.364)

  /* mRNA */
  Double("mr") default(8514.558)
  Double("mt") default(49.333)
  Double("mm") default(49.333)
  Double("mp") default(21.329)
  Double("mq") default(8811.960)

  /* ribosome-bound mRNA */
  Double("rmr") default(348.675)
  Double("rmt") default(0.0)
  Double("rmm") default(0.0)
  Double("rmp") default(0.0)
  Double("rmq") default(357.898)

  /* ribosome-bound mRNA sequestered by chloramphenicol */
  Double("zmr") default(0.0)
  Double("zmt") default(0.0)
  Double("zmm") default(0.0)
  Double("zmp") default(0.0)
  Double("zmq") default(0.0)

  // All rates should be calculated from the above.
}
