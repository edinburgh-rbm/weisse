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

import org.scalatest.{FlatSpec, Matchers}

import uk.ac.ed.inf.mois.weisse.{WeisseCell, WeisseModel}

class WeisseModelTest extends FlatSpec with Matchers {
  "an example" should "test things here" in {

    val ex = new WeisseModel

    0 should not be 1
  }
}
