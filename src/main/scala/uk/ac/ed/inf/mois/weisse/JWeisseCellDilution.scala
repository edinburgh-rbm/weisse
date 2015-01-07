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


class JWeisseCellDilution
    extends DeterministicReactionNetwork[Double, Jet[Double]]
    with Rosenbrock
    with VarCalc
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
  val lam = Double("lam")


  reactions(

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
