package uk.ac.ed.inf.mois.weisse

import uk.ac.ed.inf.mois.ode.Apache
import uk.ac.ed.inf.mois.reaction.DeterministicReactionNetwork
import uk.ac.ed.inf.mois.{Model, VarCalc, Math}
import spire.math.Jet
import spire.implicits._
import uk.ac.ed.inf.mois.implicits._

class IsolatedRepressilator
  extends DeterministicReactionNetwork[Double, Double] with Apache with VarCalc {

  val g1 = Species("g1")
  val g2 = Species("g2")
  val g3 = Species("g3")

  val mg1 = Species("mg1") default(20.0)
  val mg2 = Species("mg2")
  val mg3 = Species("mg3")

  val lambda_eff = Double("lambda_eff") default(0.022) param()
  val d_g = Double("d_g") default(0.1733) param()
  val d_mg = Double("d_mg") default(0.3466) param()
  val w_g = Double("w_g") default(500.0) param()
  val K_g = Double("K_g") default(100.0) param()
  val h = Double("h") default(2.0) param()
  val k_eff = Double("k_eff") default(0.6) param()

  val R_g1 = Double("R_g1") param()
  val R_g2 = Double("R_g2") param()
  val R_g3 = Double("R_g3") param()
  calc(R_g1) := (1.0 / (1.0 + (g3 / K_g) ** h))
  calc(R_g2) := (1.0 / (1.0 + (g1 / K_g) ** h))
  calc(R_g3) := (1.0 / (1.0 + (g2 / K_g) ** h))


  reactions(

    // Dilution
    g1 --> () at lambda_eff,
    g2 --> () at lambda_eff,
    g3 --> () at lambda_eff,
    mg1 --> () at lambda_eff,
    mg2 --> () at lambda_eff,
    mg3 --> () at lambda_eff,

    // Degradation
    g1 --> () at d_g,
    g2 --> () at d_g,
    g3 --> () at d_g,
    mg1 --> () at d_mg,
    mg2 --> () at d_mg,
    mg3 --> () at d_mg,

    // Transcription
    () --> mg1 at w_g * R_g1,
    () --> mg2 at w_g * R_g2,
    () --> mg3 at w_g * R_g3,

    // Translation
    mg1 --> mg1 + g1 at k_eff,
    mg2 --> mg2 + g2 at k_eff,
    mg3 --> mg3 + g3 at k_eff

  )
}


class IsolatedRepressilatorModel extends Model {
  val process = new IsolatedRepressilator
}
