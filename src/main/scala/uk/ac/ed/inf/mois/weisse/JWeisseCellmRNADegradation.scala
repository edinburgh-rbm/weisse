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


class JWeisseCellmRNADegradation(val dm : Double)
    extends DeterministicReactionNetwork[Double, Jet[Double]]
    with Rosenbrock
    with VarCalc
    with Math {

  /* define variables */
  /* mRNA */
  val mr = Species("mr") nonnegative()
  val mt = Species("mt") nonnegative()
  val mm = Species("mm") nonnegative()
  val mp = Species("mp") nonnegative()
  val mq = Species("mq") nonnegative()

  reactions(

    /* mRNA degradation */
    mr --> () at dm,
    mt --> () at dm,
    mm --> () at dm,
    mp --> () at dm,
    mq --> () at dm

  )
}
