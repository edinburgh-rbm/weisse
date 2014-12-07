/*
 *  Andrea Weisse' Whole Cell Model
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
name := "weisse-model"

organization := "uk.ac.ed.inf"

version := "1.99.14-SNAPSHOT"

scalaVersion := "2.11.2"

mainClass in Compile := Some("uk.ac.ed.inf.mois.MoisMain")

resolvers += "Sonatype OSS Snapshots" at "https://oss.sonatype.org/content/repositories/snapshots"

resolvers += "ucar-unidata-releases" at "https://artifacts.unidata.ucar.edu/content/repositories/unidata-releases/"

libraryDependencies += "uk.ac.ed.inf" %% "mois" % "1.99.14-SNAPSHOT"

libraryDependencies += "org.scalatest" %% "scalatest" % "2.2.0" % "test"
