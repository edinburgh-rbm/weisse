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


class JWeisseCellmRNARibosomeComplex(val kb : Double, val ku : Double)
    extends DeterministicReactionNetwork[Double, Jet[Double]]
    with Rosenbrock
    with VarCalc
    with Math {

  /* define variables */
  /* proteins */
  val r = Species("r") nonnegative()

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


  reactions(
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
    rmq --> r + mq at ku

  )
}
