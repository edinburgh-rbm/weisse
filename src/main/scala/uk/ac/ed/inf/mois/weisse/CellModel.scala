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

import uk.ac.ed.inf.mois.{Model, Process}
import uk.ac.ed.inf.mois.{DeterministicReactionNetwork, VarCalc, Math}

import org.apache.commons.math3.ode.nonstiff.AdamsMoultonIntegrator

class WeisseCell(
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
    extends DeterministicReactionNetwork("AndreasCell")
       with VarCalc
       with Math {

  override def integrator = {
    val nSteps = 4
    val minStep = 1e-8
    val maxStep = 100
    val absoluteTolerance = 1e-10
    val relativeTolerance = 1e-10
    new AdamsMoultonIntegrator(nSteps, minStep, maxStep,
      absoluteTolerance, relativeTolerance)
  }

  import MultisetConversion._

  /* define variables */
  /* ATP and internal nutrient */
  val a = Species("a")
  val si = Species("si")

  /* proteins */
  val r = Species("r")
  val et = Species("et")
  val em = Species("em")
  val p = Species("p")
  val q = Species("q")

  /* mRNA */
  val mr = Species("mr")
  val mt = Species("mt")
  val mm = Species("mm")
  val mp = Species("mp")
  val mq = Species("mq")

  /* ribosome-bound mRNA */
  val rmr = Species("rmr")
  val rmt = Species("rmt")
  val rmm = Species("rmm")
  val rmp = Species("rmp")
  val rmq = Species("rmq")

  /* ribosome-bound mRNA sequestered by chloramphenicol */
  val zmr = Species("zmr")
  val zmt = Species("zmt")
  val zmm = Species("zmm")
  val zmp = Species("zmp")
  val zmq = Species("zmq")

  /* some rates */
  val Kgamma = Double("Kgamma")
  val gamma = Double("gamma")
  val ttrate = Double("ttrate")
  val lam = Double("lam")
  val fr = Double("fr")
  val nucat = Double("nucat")
  val f = Double("f")
  val b = Double("b")

  calc(Kgamma) := gmax/Kp
  calc(gamma) := gmax*a/(Kgamma + a)
  calc(ttrate) := gamma*(rmq + rmr + rmp + rmt + rmm)
  calc(lam) := ttrate/M
  calc(fr) := nr*(r + rmr + rmp + rmt + rmm + rmq + zmr + zmp + zmt + zmm + zmq) /
    ( nr*(r + rmr + rmp + rmt + rmm + rmq + zmr + zmp + zmt + zmm + zmq) + nx*(p + q + et + em) )
  calc(nucat) := em*vm*si/(Km + si)
  calc(f) := cl*k_cm
  calc(b) := 0.0

  reactions(
    /* nutrient import */
    () -> si `at!` et*vt*s0/(Kt + s0),

    /* nutrient metabolism */
    si -> () `at!` nucat,
    () -> a `at!` ns*nucat,
    a -> () `at!` ttrate,

    /* transcription */
    () -> mr `at!` wr*a/(thetar + a),
    () -> mt `at!` we*a/(thetax + a),
    () -> mm `at!` we*a/(thetax + a),
    () -> mp `at!` wp*a/(thetax + a),
    () -> mq `at!` wq*a/(thetax + a)/(1 + pow((q/Kq), nq)),

    /* translation */
    r + mr -> rmr at kb,
    r + mt -> rmt at kb,
    r + mm -> rmm at kb,
    r + mp -> rmp at kb,
    r + mq -> rmq at kb,

    rmr -> r + mr at ku,
    rmt -> r + mt at ku,
    rmm -> r + mm at ku,
    rmp -> r + mp at ku,
    rmq -> r + mq at ku,

    rmr -> r + r + mr `at!` gamma/nr,
    rmt -> r + et + mt `at!` gamma/nx,
    rmm -> r + em + mm `at!` gamma/nx,
    rmp -> r + p + mp `at!` gamma/nx,
    rmq -> r + q + mq `at!` gamma/nx,

    /* mRNA degradation */
    mr -> () at dm,
    mt -> () at dm,
    mm -> () at dm,
    mp -> () at dm,
    mq -> () at dm,

    /* chloramphenicol */
    rmr -> zmr at f,
    rmt -> zmt at f,
    rmm -> zmm at f,
    rmp -> zmp at f,
    rmq -> zmq at f,

    zmr -> rmr at b,
    zmt -> rmt at b,
    zmm -> rmm at b,
    zmp -> rmp at b,
    zmq -> rmq at b,

    /* dilution */
    mr -> () at lam,
    mt -> () at lam,
    mm -> () at lam,
    mp -> () at lam,
    mq -> () at lam,
    rmr -> () at lam,
    rmt -> () at lam,
    rmm -> () at lam,
    rmp -> () at lam,
    rmq -> () at lam,
    r -> () at lam,
    et -> () at lam,
    em -> () at lam,
    p -> () at lam,
    q -> () at lam,
    a -> () at lam,
    si -> () at lam,
    zmr -> () at lam,
    zmt -> () at lam,
    zmm -> () at lam,
    zmp -> () at lam,
    zmq -> () at lam
  )
}

class WeisseModel extends Model {

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

  val process = new WeisseCell(
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

  import process._

  // -- Variable constraints --

  /* ATP and internal nutrient */
  a must (_ >= 0.0)
  si must (_ >= 0.0)

  /* proteins */
  r must (_ >= 0.0)
  et must (_ >= 0.0)
  em must (_ >= 0.0)
  p must (_ >= 0.0)
  q must (_ >= 0.0)

  /* mRNA */
  mr must (_ >= 0.0)
  mt must (_ >= 0.0)
  mm must (_ >= 0.0)
  mp must (_ >= 0.0)
  mq must (_ >= 0.0)

  /* ribosome-bound mRNA */
  rmr must (_ >= 0.0)
  rmt must (_ >= 0.0)
  rmm must (_ >= 0.0)
  rmp must (_ >= 0.0)
  rmq must (_ >= 0.0)

  /* ribosome-bound mRNA sequestered by chloramphenicol */
  zmr must (_ >= 0.0)
  zmt must (_ >= 0.0)
  zmm must (_ >= 0.0)
  zmp must (_ >= 0.0)
  zmq must (_ >= 0.0)

  /* some rates */
  Kgamma must (_ >= 0.0)
  gamma must (_ >= 0.0)
  ttrate must (_ >= 0.0)
  lam must (_ >= 0.0)
  fr must (_ >= 0.0)
  nucat must (_ >= 0.0)
  f must (_ >= 0.0)
  b must (_ >= 0.0)

  // -- Initial values --

  a := 1.0
  r := 1.0

  // Everything else is just zero
}
