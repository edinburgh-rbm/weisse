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


class JWeisseCellChloramphenicol
    extends DeterministicReactionNetwork
       with VarCalc
       with Math {

  /* define variables */
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
  val f = Double("f")
  val b = Double("b")


  reactions(

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
    zmq --> rmq at b

  )
}