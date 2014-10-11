 /* 
  *  MOIS Examples: Mammalian Circadian Clock (Forger-Peskin 2002)
  *  Copyright (C) 2014 University of Edinburgh School of Informatics
  * 
  *  This program is free software: you can redistribute it and/or modify
  *  it under the terms of the GNU General Public License as published by
  *  the Free Software Foundation, either version 3 of the License, or
  *  (at your option) any later version.
  * 
  *  This program is distributed in the hope that it will be useful,
  *  but WITHOUT ANY WARRANTY; without even the implied warranty of
  *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  *  GNU General Public License for more details.
  * 
  *  You should have received a copy of the GNU General Public License
  *  along with this program. If not, see <http://www.gnu.org/licenses/>.
  */
package uk.ac.ed.inf.mois.clocks.test

import org.scalatest.{FlatSpec, Matchers}

import uk.ac.ed.inf.mois.clocks.{MammalianCircadianClock,MammalianCircadianClockModel}

class MammalianCircadianClockModelTest extends FlatSpec with Matchers {
  "an example" should "test things here" in {

    val ex = new MammalianCircadianClockModel

  /* Both quantities below should remain constant during simulation:
   *
   *  C = Ct - (PoC + PtC + PopC + PtpC + PoppC + PtppC + PopCRo + PopCRt + PtpCRo + PtpCRt + PoppCRo + PoppCRt + PtppCRo + PtppCRt + PonpCn + PtnpCn + PonppCn
+ PtnppCn + PonpCnRon + PonpCnRtn + PtnpCnRon + PtnpCnRtn + PonppCnRon + PonppCnRtn + PtnppCnRon + PtnppCnRtn + Cn)
   *
   *  Rn = (Ron + PonpRon + PonppRon + PonpCnRon + PonppCnRon + PtnpRon + PtnppRon + PtnpCnRon + PtnppCnRon + Rtn + PonpRtn + PonppRtn + PonpCnRtn + PonppCnRtn + PtnpRtn + PtnppRtn + PtnpCnRtn + PtnppCnRtn
   *
   */

  }
}
