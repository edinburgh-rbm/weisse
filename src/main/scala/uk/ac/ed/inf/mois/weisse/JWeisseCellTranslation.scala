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


class JWeisseCellTranslation(val nr : Double, val nx : Double)
    extends DeterministicReactionNetwork
       with VarCalc
       with Math {

  /* define variables */
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

  /* some rates */
  val gamma = Double("gamma")


  reactions(

    /* translation */
    rmr --> r + r + mr `at!` gamma/nr,
    rmt --> r + et + mt `at!` gamma/nx,
    rmm --> r + em + mm `at!` gamma/nx,
    rmp --> r + p + mp `at!` gamma/nx,
    rmq --> r + q + mq `at!` gamma/nx

  )
}