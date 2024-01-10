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

package at.tugraz.genome.lda.swing;

import java.awt.Color;
import java.awt.Font;
import java.awt.event.ActionListener;

import javax.swing.JLabel;
import javax.swing.JPanel;

import at.tugraz.genome.lda.Settings;
import at.tugraz.genome.lda.TooltipTexts;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class ExportPanel extends JPanel
{
  
  private static final long serialVersionUID = 1684790632349952418L;
  
  public static final String EXPORT_SVG = "exportSVG";
  
  public static final String EXPORT_PNG = "exportPNG";
  
  public static final String EXPORT_EXCEL = "exportExcel";
  
  public static final String EXPORT_TEXT = "exportText";
  
  public static final String EXPORT_CHROMS = "exportChroms";
  
  public static final String EXPORT_MZTAB = "exportMztab";
  
  public static final String EXPORT_RDB = "exportRdb";

  public static final String EXPORT_MAF = "exportMAF";
  
  public static final String EXPORT_SUMMARY = "exportSummary";

  
  private ActionListener parent_;
  
  public ExportPanel(Color background, Color font, ActionListener parent){
    this(background,font,parent,false,false,false);
  }
  
  public ExportPanel(Color background, Color font, ActionListener parent, boolean chroms, boolean mzTab, boolean pictureOnly){
    this.parent_ = parent;
    if (background != null)
      this.setBackground(background);
    Font textFont = new Font("Helvetica",Font.BOLD, 9);

    JLabel label = null;
    if (!pictureOnly){
      label = new JLabel("Export: ");
      label.setForeground(font);
      label.setFont(textFont);
      label.setToolTipText(TooltipTexts.EXPORT_GENERAL);
      this.add(label);
    }
    ExportButton png = new ExportButton("PNG",EXPORT_PNG,font,background,parent_);
    png.setToolTipText(TooltipTexts.EXPORT_PNG);
    this.add(png);
    label = new JLabel(" | ");
    label.setForeground(font);
    label.setFont(textFont);
    this.add(label); 
    ExportButton svg = new ExportButton("SVG",EXPORT_SVG,font,background,parent_);
    svg.setToolTipText(TooltipTexts.EXPORT_SVG);
    this.add(svg);
//    if (pictureOnly){
//      label = new JLabel(" | ");
//      label.setForeground(font);
//      label.setFont(textFont);
//      this.add(label); 
//      ExportButton mgf = new ExportButton("mgf","mgf",font,background,parent_);
//      //mgf.setToolTipText(TooltipTexts.EXPORT_SVG);
//      this.add(mgf);     
//    }
    if (pictureOnly) return;
    label = new JLabel(" | ");
    label.setForeground(font);
    label.setFont(textFont);
    this.add(label); 
    ExportButton excel = new ExportButton("Excel",EXPORT_EXCEL,font,background,parent_);
    excel.setToolTipText(TooltipTexts.EXPORT_EXCEL);
    this.add(excel);
    if (mzTab){
      label = new JLabel(" | ");
      label.setForeground(font);
      label.setFont(textFont);
      this.add(label); 
      ExportButton mztab = new ExportButton("mzTab",EXPORT_MZTAB,font,background,parent_);
      mztab.setToolTipText(TooltipTexts.EXPORT_MZTAB);
      this.add(mztab);
      label = new JLabel(" | ");
      label.setForeground(font);
      label.setFont(textFont);
      this.add(label); 
      ExportButton rdb = new ExportButton("rdb",EXPORT_RDB,font,background,parent_);
      rdb.setToolTipText(TooltipTexts.EXPORT_RDB);
      this.add(rdb);
    }
    label = new JLabel(" | ");
    label.setForeground(font);
    label.setFont(textFont);
    this.add(label); 
    ExportButton text = new ExportButton("Text",EXPORT_TEXT,font,background,parent_);
    text.setToolTipText(TooltipTexts.EXPORT_TXT);
    this.add(text);
    if (chroms){
      label = new JLabel(" | ");
      label.setForeground(font);
      label.setFont(textFont);
      this.add(label); 
      ExportButton chrom = new ExportButton("Chroms",EXPORT_CHROMS,font,background,parent_);
      chrom.setToolTipText(TooltipTexts.EXPORT_CHROMS);
      this.add(chrom);      
    }
    if (Settings.SHOW_OMEGA_TOOLS) //TODO: for SILDA analysis, remove when done
    {
      label = new JLabel(" | ");
      label.setForeground(font);
      label.setFont(textFont);
      this.add(label); 
      ExportButton summary = new ExportButton("Summary",EXPORT_SUMMARY,font,background,parent_);
      this.add(summary);      
    }
    
    

    //// this has to be removed in the default LDA version; START!
//    label = new JLabel(" | ");
//    label.setForeground(font);
//    label.setFont(textFont);
//    this.add(label); 
//    ExportButton maf = new ExportButton("MAF",EXPORT_MAF,font,background,parent_);
//    maf.setToolTipText(TooltipTexts.EXPORT_MAF);
//    this.add(maf);
    ////this has to be removed in the default LDA version; END!
  }

}
