/*
 *  Tests for Andrea Weisse' cell model
 *  Copyright (C) 2014 Weisse et al
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
package uk.ac.ed.inf.mois.weisse.test

import uk.ac.ed.inf.mois.implicits._

import scala.language.reflectiveCalls
import org.scalatest.{FlatSpec, Matchers}
import org.scalactic.TolerantNumerics

import uk.ac.ed.inf.mois.{Model, ProcessGroup, State}
import uk.ac.ed.inf.mois.ode.{ODE, Apache}

import uk.ac.ed.inf.mois.weisse._

class WeisseModelSteadyTest extends FlatSpec with Matchers {
  "Andrea's equalities" should "hold true" in {

  // Use approximate equality in `should equal`
    val precision = 1e-3
    implicit val doubleEquality =
     TolerantNumerics.tolerantDoubleEquality(precision)

    val model = new WeisseModelSteady
    import model._
    model.init(0)

    val s = 1e-4
    for (i <- 0.0 until 1.0 by s) {
      model.process.step(i, i+s)

      (em.value * vm / (Km + si.value)) should equal (et.value * vt / (Kt + s0))
    }
  }
}
