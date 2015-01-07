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
import uk.ac.ed.inf.mois.ode.{Apache, Rosenbrock}
import uk.ac.ed.inf.mois.reaction.DeterministicReactionNetwork
import spire.math.Jet
import spire.implicits._
import uk.ac.ed.inf.mois.implicits._


class JWeisseCellTranscription(val thetar: Double, val thetax: Double, val wr : Double, val wq : Double, val we : Double, val wp : Double, val Kq : Double, val nq : Double)
    extends DeterministicReactionNetwork[Double, Double]
    with Apache
    with VarCalc
    with Math {

  /* define variables */
  /* ATP and internal nutrient */
  val a = Species("a") nonnegative()

  /* proteins */
  val q = Species("q") nonnegative()

  /* mRNA */
  val mr = Species("mr") nonnegative()
  val mt = Species("mt") nonnegative()
  val mm = Species("mm") nonnegative()
  val mp = Species("mp") nonnegative()
  val mq = Species("mq") nonnegative()


  reactions(

    /* transcription */
    () --> mr `at!` wr*a/(thetar + a),
    () --> mt `at!` we*a/(thetax + a),
    () --> mm `at!` we*a/(thetax + a),
    () --> mp `at!` wp*a/(thetax + a),
    () --> mq `at!` wq*a/(thetax + a)/(1 + pow((q/Kq), nq))

  )
}
