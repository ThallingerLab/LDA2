/* 
 * This file is part of Lipid Data Analyzer
 * Lipid Data Analyzer - Automated annotation of lipid species and their molecular structures in high-throughput data from tandem mass spectrometry
 * Copyright (c) 2017 Juergen Hartler, Andreas Ziegl, Gerhard G. Thallinger, Leonida M. Lamp 
 * DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS FILE HEADER. 
 *  
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * by the Free Software Foundation, either version 3 of the License, or 
 * (at your option) any later version.
 *  
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details. 
 *  
 * You should have received a copy of the GNU General Public License 
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 *
 * Please contact lda@genome.tugraz.at if you need additional information or 
 * have any questions.
 */ 

package at.tugraz.genome.lda.xml;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class XMLConstants
{
  public static final String WHOLE_DTD_NAME = "WholeSet";
  
  public static final String WHOLE_ABS_SET_ROOT = "abs-settings-whole";
  public static final String WHOLE_ABS_SET_EXPS_SET = "exp-settings";
  public static final String WHOLE_ABS_SET_EXP_SET = "exp-setting";
  public static final String WHOLE_ABS_SET_EXP_SET_ATTR_NAME = "name";
  
  public static final String WHOLE_ABS_NAME_PROBE_VOL="probe-volume";
  public static final String WHOLE_ABS_NAME_END_VOL="end-volume";
  public static final String WHOLE_ABS_NAME_SAMPLE_WEIGHT="sample-weight";

  public static final String WHOLE_ABS_NAME_PROTEIN_CONC="protein-concentration";
  public static final String WHOLE_ABS_NAME_LIPID_CONC="neutral-lipid-concentration";
  public static final String WHOLE_ABS_NAME_DILUTION="dilution-factor";
  public static final String WHOLE_ABS_NAME_ALL_EXPS_SAME="all-experiments-same";
  public static final String WHOLE_ABS_NAME_ALL_STANDS_SAME="all-standards-same";
  public static final String WHOLE_ABS_NAME_EXT_VOLUME="ext-volume";
  public static final String WHOLE_ABS_NAME_EXT_CONC="ext-concentration";
  public static final String WHOLE_ABS_NAME_INT_VOLUME="int-volume";
  public static final String WHOLE_ABS_NAME_INT_CONC="int-concentration";
  
  public static final String WHOLE_ABS_SET_EXPS_STAND_SET = "exp-stand-settings";
  public static final String WHOLE_ABS_SET_EXP_STAND_SET = "exp-stand-setting";
  
  public static final String ES_SETTINGS_VALUE_GENERAL = "general";
  public static final String IS_SETTINGS_VALUE_GENERAL = "general";
  
  public static final String WHOLE_ABS_SET_CLASSES_SET = "class-settings";
  public static final String WHOLE_ABS_SET_CLASS_SET = "class-setting";
  public static final String WHOLE_ABS_SET_CLASS_NAME = "lipid-class";
  
  public static final String TAG_SETTING = "setting";
  public static final String TAG_ES_SETTING = "es-setting";
  public static final String TAG_IS_SETTING = "is-setting";
  
  public static final String TAG_SETTING_ATTR_MAGNITUDE = "magnitude";
  
  public static final String ATTR_NAME = "name";
  public static final String ATTR_VALUE = "value";
  public static final String ATTR_CLASS = "class";
  
  
  public static final String CUTOFF_DTD_NAME = "CutoffSet";
  
  public static final String CUTOFF_SET_ROOT = "cutoff-settings";
  public static final String CUTOFF_ISOTOPE = "isotope";
  public static final String CUTOFF_CUTOFF = "cutoff";
  
}
