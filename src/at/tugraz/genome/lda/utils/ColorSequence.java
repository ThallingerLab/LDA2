/* 
 * This file is part of Lipid Data Analyzer
 * Lipid Data Analyzer - Automated annotation of lipid species and their molecular structures in high-throughput data from tandem mass spectrometry
 * Copyright (c) 2017 Juergen Hartler, Andreas Ziegl, Gerhard G. Thallinger 
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

package at.tugraz.genome.lda.utils;

import java.awt.Color;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class ColorSequence
{
  public static final Color VERY_DARK_RED = new Color(0x80, 0x00, 0x00);
  public static final Color DARK_RED = new Color(0xc0, 0x00, 0x00);
  public static final Color LIGHT_RED = new Color(0xFF, 0x40, 0x40);
  public static final Color VERY_LIGHT_RED = new Color(0xFF, 0x80, 0x80);
  public static final Color VERY_DARK_YELLOW = new Color(0x80, 0x80, 0x00);
  public static final Color DARK_YELLOW = new Color(0xC0, 0xC0, 0x00);
  public static final Color LIGHT_YELLOW = new Color(0xFF, 0xFF, 0x40);
  public static final Color VERY_LIGHT_YELLOW = new Color(0xFF, 0xFF, 0x80);
  public static final Color VERY_DARK_GREEN = new Color(0x00, 0x80, 0x00);
  public static final Color DARK_GREEN = new Color(0x00, 0xC0, 0x00);
  public static final Color LIGHT_GREEN = new Color(0x40, 0xFF, 0x40);
  public static final Color VERY_LIGHT_GREEN = new Color(0x80, 0xFF, 0x80);
  public static final Color VERY_DARK_CYAN = new Color(0x00, 0x80, 0x80);
  public static final Color DARK_CYAN = new Color(0x00, 0xC0, 0xC0);
  public static final Color LIGHT_CYAN = new Color(0x40, 0xFF, 0xFF);
  public static final Color VERY_LIGHT_CYAN = new Color(0x80, 0xFF, 0xFF);
  public static final Color VERY_DARK_BLUE = new Color(0x00, 0x00, 0x80);
  public static final Color DARK_BLUE = new Color(0x00, 0x00, 0xC0);
  public static final Color LIGHT_BLUE = new Color(0x40, 0x40, 0xFF);
  public static final Color VERY_LIGHT_BLUE = new Color(0x80, 0x80, 0xFF);
  public static final Color VERY_DARK_MAGENTA = new Color(0x80, 0x00, 0x80);
  public static final Color DARK_MAGENTA = new Color(0xC0, 0x00, 0xC0);
  public static final Color LIGHT_MAGENTA = new Color(0xFF, 0x40, 0xFF);
  public static final Color VERY_LIGHT_MAGENTA = new Color(0xFF, 0x80, 0xFF);

  
  public static Color[] getDefaultColors(){
    return new Color[]{new Color(0xFF, 0x55, 0x55),
        new Color(0x55, 0x55, 0xFF),
        new Color(0x55, 0xFF, 0x55),
        new Color(0xFF, 0xFF, 0x55),
        new Color(0xFF, 0x55, 0xFF),
        new Color(0x55, 0xFF, 0xFF),
        Color.pink,
        Color.gray,
        ColorSequence.LIGHT_RED,
        ColorSequence.LIGHT_BLUE,
        ColorSequence.LIGHT_GREEN,
        ColorSequence.LIGHT_YELLOW,
        ColorSequence.LIGHT_MAGENTA,
        ColorSequence.LIGHT_CYAN,
        Color.lightGray,
        ColorSequence.VERY_DARK_RED,
        ColorSequence.VERY_DARK_BLUE,
        ColorSequence.VERY_DARK_GREEN,
        ColorSequence.VERY_DARK_YELLOW,
        ColorSequence.VERY_DARK_MAGENTA,
        ColorSequence.VERY_DARK_CYAN,
        ColorSequence.VERY_LIGHT_RED,
        ColorSequence.VERY_LIGHT_BLUE,
        ColorSequence.VERY_LIGHT_GREEN,
        ColorSequence.VERY_LIGHT_YELLOW,
        ColorSequence.VERY_LIGHT_MAGENTA,
        ColorSequence.VERY_LIGHT_CYAN,
        ColorSequence.DARK_RED,
        ColorSequence.DARK_BLUE,
        ColorSequence.DARK_GREEN,
        ColorSequence.DARK_YELLOW,
        ColorSequence.DARK_MAGENTA,
        ColorSequence.DARK_CYAN,
        Color.darkGray};
  }
  
}
