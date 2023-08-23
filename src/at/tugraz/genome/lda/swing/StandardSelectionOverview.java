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

import java.awt.BorderLayout;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionListener;
import java.util.Hashtable;

import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JPanel;

import at.tugraz.genome.lda.TooltipTexts;
import at.tugraz.genome.lda.vos.ResultCompVO;
import at.tugraz.genome.lda.vos.ResultDisplaySettingsVO;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class StandardSelectionOverview extends JDialog
{

  private static final long serialVersionUID = -6125157095806596556L;

  private Hashtable<String,ResultDisplaySettings> displayNames_;
  private ActionListener parent_;
  private Hashtable<String,JLabel> labelTexts_;
  private Hashtable<String,Hashtable<String,Integer>> isLookup_;
  private Hashtable<String,Hashtable<String,Integer>> esLookup_;
  
  private JButton buttonOK_;
  private JButton buttonCancel_;
  private JPanel informationPanel_;
  
  public StandardSelectionOverview(Hashtable<String,ResultDisplaySettings> displayNames, ActionListener parent,
      Hashtable<String,Hashtable<String,Integer>> isLookup, Hashtable<String,Hashtable<String,Integer>> esLookup){
    displayNames_ = displayNames;
    parent_ = parent;
    isLookup_ = isLookup;
    esLookup_ = esLookup;
    setLayout(new BorderLayout());
    setLocation(380,240);
    initOverview();
    initButtonPanel();
    this.setVisible(false);
    setDefaultCloseOperation(DISPOSE_ON_CLOSE);
    pack(); 
  }
  
  private void initOverview(){
    JLabel topPanel = new JLabel("The following standards are used for comparison of the classes");
    add(topPanel,BorderLayout.NORTH);
    informationPanel_ = new JPanel();
    informationPanel_.setLayout(new GridBagLayout());
    labelTexts_ = new Hashtable<String,JLabel>();
    int count = 0;
    for (String className : displayNames_.keySet()){
      JLabel classLabel = new JLabel(className+": ");
      informationPanel_.add(classLabel,new GridBagConstraints(0, count, 1, 1, 0.0, 0.0
          ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 2, 0, 2), 0, 0));
      JLabel infoLabel = new JLabel("");
      informationPanel_.add(infoLabel,new GridBagConstraints(1, count, 1, 1, 0.0, 0.0
          ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 2, 0, 2), 0, 0));
      labelTexts_.put(className, infoLabel);
      count++;
    }
    add(informationPanel_,BorderLayout.CENTER);
  }
  
  private void initButtonPanel(){
    JPanel buttonPanel = new JPanel();
    buttonPanel.setLayout(new GridBagLayout());
    buttonOK_ = new JButton("OK");
    buttonOK_.addActionListener(parent_);
    buttonOK_.setActionCommand("acceptSelectedStandards");
    buttonOK_.setToolTipText(TooltipTexts.ACCEPT_GENERAL);
    buttonPanel.add(buttonOK_,new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 2, 0, 2), 0, 0));
    buttonCancel_ = new JButton("Cancel");
    buttonCancel_.addActionListener(parent_);
    buttonCancel_.setActionCommand("cancelSelectedStandards");
    buttonCancel_.setToolTipText(TooltipTexts.CANCEL_GENERAL);
    buttonPanel.add(buttonCancel_,new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 2, 0, 2), 0, 0));
    add(buttonPanel,BorderLayout.SOUTH);
  }
  
  public void refreshInformationSections(){
    for (String className : displayNames_.keySet()){
      ResultDisplaySettingsVO settings = displayNames_.get(className).getSettingsVO();
      JLabel infoLabel = labelTexts_.get(className);
      String labelText = "";
      if (settings.getISStandMethod()==ResultCompVO.NO_STANDARD_CORRECTION && settings.getESStandMethod()==ResultCompVO.NO_STANDARD_CORRECTION){
        labelText = "no standard available - class is neglected!";
      }else{
        if (settings.getISStandMethod()!=ResultCompVO.NO_STANDARD_CORRECTION){
          labelText += "IS: ";
          if (settings.getISStandMethod()==ResultCompVO.STANDARD_CORRECTION_INTERNAL){
            labelText += "most reliable";
          } else if (settings.getISStandMethod()==ResultCompVO.STANDARD_CORRECTION_MEDIAN) {
            labelText += "median";
          } else {
            Hashtable<String,Integer> isLook = isLookup_.get(className);
            for (String name : isLook.keySet()){
              if (settings.getISStandMethod()==isLook.get(name)){
                labelText += name;
                break;
              }
            }
          }
          if (settings.getESStandMethod()!=ResultCompVO.NO_STANDARD_CORRECTION){
            labelText += " -> ";
          }
        }
        if (settings.getESStandMethod()!=ResultCompVO.NO_STANDARD_CORRECTION){
          labelText += "ES: ";
          if (settings.getESStandMethod()==ResultCompVO.STANDARD_CORRECTION_INTERNAL){
            labelText += "most reliable";
          } else if (settings.getESStandMethod()==ResultCompVO.STANDARD_CORRECTION_MEDIAN) {
            labelText += "median";
          } else {
            Hashtable<String,Integer> esLook = esLookup_.get(className);
            for (String name : esLook.keySet()){
              if (settings.getESStandMethod()==esLook.get(name)){
                labelText += name;
                break;
              }
            }
          }
        }
      }
      
      infoLabel.setText(labelText);
    }
    invalidate();
    doLayout();
    this.setSize(this.getPreferredSize());
    repaint();
  }
  
}

