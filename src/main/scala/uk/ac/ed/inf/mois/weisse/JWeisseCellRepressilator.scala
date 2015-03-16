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
import uk.ac.ed.inf.mois.ode.Rosenbrock
import uk.ac.ed.inf.mois.reaction.DeterministicReactionNetwork
import spire.math.Jet
import spire.implicits._
import uk.ac.ed.inf.mois.implicits._


class JWeisseCellRepressilator(val kb: Double, val ku: Double, val thetax: Double)
    extends DeterministicReactionNetwork[Double, Jet[Double]]
    with Rosenbrock
    with VarCalc
    with Math {

  /* define variables */
  /* ATP and internal nutrient */
  val a = Species("a") nonnegative()

  /* proteins */
  val r = Species("r") nonnegative()

  /* some rates */
  val gamma = Double("gamma")
  val lam = Double("lam")


  /* Repressilator species */
  val g1 = Species("g1")
  val g2 = Species("g2")
  val g3 = Species("g3")

  val mg1 = Species("mg1") default(20.0)
  val mg2 = Species("mg2")
  val mg3 = Species("mg3")

  val cg1 = Species("cg1")
  val cg2 = Species("cg2")
  val cg3 = Species("cg3")

  /* Repressilator rates */

  val d_g = Double("d_g") default(0.1733) param()
  val d_mg = Double("d_mg") default(0.3466) param()
  val w_g = Double("w_g") default(500.0) param()
  val K_g = Double("K_g") default(100.0) param()
  val h = Double("h") default(2.0) param()
  val n_g = Double("n_g") default(300.0) param()

  import spire.math.Jet
  implicit def jetToDouble(j: Jet[Double]): Double = j.real

  val R_g1 = Double("R_g1") param()
  val R_g2 = Double("R_g2") param()
  val R_g3 = Double("R_g3") param()
  calc(R_g1) := (1.0 / (1.0 + (jetToDouble(g3) / K_g) ** h))
  calc(R_g2) := (1.0 / (1.0 + (jetToDouble(g1) / K_g) ** h))
  calc(R_g3) := (1.0 / (1.0 + (jetToDouble(g2) / K_g) ** h))

  val w_g1 = Double("w_g1") param()
  val w_g2 = Double("w_g2") param()
  val w_g3 = Double("w_g3") param()
  calc(w_g1) := (w_g * ((a / (thetax + a)) * R_g1))
  calc(w_g2) := (w_g * ((a / (thetax + a)) * R_g2))
  calc(w_g3) := (w_g * ((a / (thetax + a)) * R_g3))

  val rttrate = Double("rttrate") param()
  calc(rttrate) := gamma / n_g


  reactions(

    // Dilution
    g1 --> () at lam,
    g2 --> () at lam,
    g3 --> () at lam,
    mg1 --> () at lam,
    mg2 --> () at lam,
    mg3 --> () at lam,
    cg1 --> () at lam,
    cg2 --> () at lam,
    cg3 --> () at lam,

    // Degradation
    g1 --> () at d_g,
    g2 --> () at d_g,
    g3 --> () at d_g,
    mg1 --> () at d_mg,
    mg2 --> () at d_mg,
    mg3 --> () at d_mg,

    // Transcription
    () --> mg1 at w_g1,
    () --> mg2 at w_g2,
    () --> mg3 at w_g3,

    // Ribosome Binding
    r + mg1 --> cg1 at kb,
    r + mg2 --> cg2 at kb,
    r + mg3 --> cg3 at kb,

    cg1 --> r + mg1 at ku,
    cg2 --> r + mg2 at ku,
    cg3 --> r + mg3 at ku,

    // Translation
    cg1 --> r + mg1 + g1 at rttrate,
    cg2 --> r + mg2 + g2 at rttrate,
    cg3 --> r + mg3 + g3 at rttrate

  )
}
