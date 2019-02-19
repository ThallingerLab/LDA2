/* 
 * This file is part of Lipid Data Analyzer
 * Lipid Data Analyzer - Automated annotation of lipid species and their molecular structures in high-throughput data from tandem mass spectrometry
 * Copyright (c) 2018 Juergen Hartler 
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
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Vector;

import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JPanel;

import at.tugraz.genome.lda.TooltipTexts;
import at.tugraz.genome.lda.exception.LipidCombinameEncodingException;
import at.tugraz.genome.lda.msn.LipidomicsMSnSet;
import at.tugraz.genome.lda.msn.MSnAnalyzer;
import at.tugraz.genome.lda.quantification.LipidParameterSet;
import at.tugraz.genome.lda.vos.DoubleStringVO;
import at.tugraz.genome.maspectras.utils.Calculator;
import at.tugraz.genome.voutils.GeneralComparator;

/**
 * displays the recalculated MSn identification
 * 
 * @author Juergen Hartler
 *
 */
public class RecalculateMSnDialog extends JDialog implements ActionListener
{

  private static final long serialVersionUID = -7564408227460237486L;

  /** the detected MSn result*/
  private LipidParameterSet result_;
  /** the parent action listener*/
  private ActionListener parent_;
  /** is an accept of the result allowed*/
  private boolean acceptAllowed_;
  /** the analyte class*/
  private String className_;

  /**
   * Constructor for a dialog field showing the new MSn assignment
   * @param className the name of the lipid class
   * @param analyzer the MSnAnalyzer that calculated the new assignment
   * @param parent the parent listener that is informed when a button is pressed
   */
  public RecalculateMSnDialog(String className, MSnAnalyzer analyzer, ActionListener parent){
    acceptAllowed_ = true;
    this.className_ = className;
    this.parent_ = parent;
    this.setLayout(new BorderLayout());
    setLocation(380,240);
    result_ = analyzer.getResult();
    initMainPanel(analyzer);
    initButtonPanel();
    pack(); 
    setVisible(true);
  }
  
  /**
   * initializes the main panel that contains the input field for the retention time
   * @param analyzer the MSnAnalyzer that calculated the new assignment
   */
  @SuppressWarnings("unchecked")
  private void initMainPanel(MSnAnalyzer analyzer){
    int status = analyzer.checkStatus();
    JPanel centerPanel = new JPanel();
    centerPanel.setLayout(new GridBagLayout());
    JLabel label;
    if (status==LipidomicsMSnSet.NO_MSN_PRESENT){
      label = new JLabel("There are no MSn spectra present!");
      centerPanel.add(label,new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
          ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(1, 1, 1, 1), 0, 0));
      label = new JLabel("Previous MSn detections will be deleted when you click accept");
      centerPanel.add(label,new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0
          ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(1, 1, 1, 1), 0, 0));
    } else if (status==LipidomicsMSnSet.DISCARD_HIT){
      label = new JLabel("The MSn recommendation is to reject this hit!");
      centerPanel.add(label,new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
          ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(1, 1, 1, 1), 0, 0));
      label = new JLabel("Please delete it yourself!");
      centerPanel.add(label,new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0
          ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(1, 1, 1, 1), 0, 0));
      acceptAllowed_ = false;
    } else if (status==LipidomicsMSnSet.HEAD_GROUP_DETECTED){
      label = new JLabel("The species could be verified as "+className_+" "+result_.getNameStringWithoutRt()+"!");
      centerPanel.add(label,new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
          ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(1, 1, 1, 1), 0, 0));
      label = new JLabel("There were no molecular species detectable");
      centerPanel.add(label,new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0
          ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(1, 1, 1, 1), 0, 0));
      label = new JLabel("Please click accept when you want to store this recommendation");
      centerPanel.add(label,new GridBagConstraints(0, 2, 1, 1, 0.0, 0.0
          ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(1, 1, 1, 1), 0, 0));
    } else if (status>LipidomicsMSnSet.HEAD_GROUP_DETECTED){
      label = new JLabel("The molecular species were identified:");
      centerPanel.add(label,new GridBagConstraints(0, 0, 3, 1, 0.0, 0.0
          ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(1, 1, 1, 1), 0, 0));
      float area = result_.getArea();
      LipidomicsMSnSet msn = (LipidomicsMSnSet)result_;
      List<DoubleStringVO> nameAreaVO = new ArrayList<DoubleStringVO>();
      Vector<Object> detected = null;
      try {detected = msn.getMSnIdentificationNames();
      }catch (LipidCombinameEncodingException lcx) {
        detected = new Vector<Object>();
        lcx.printStackTrace();
      }
      for (Object names : detected){
        String name = "";
        double relArea = 0d;
        if (names instanceof Vector){
          Vector<String> nameSuggestions = (Vector<String>)names;
          relArea =msn.getRelativeIntensity(nameSuggestions.get(0));
          for (String nameSuggestion : nameSuggestions){
            if (name.length()>0) name+=";";
            name += nameSuggestion;
          }
        } else if (names instanceof String){
          name = (String)names;
          relArea =msn.getRelativeIntensity(name);
        }
        nameAreaVO.add(new DoubleStringVO(name,relArea));
      }
      Collections.sort(nameAreaVO,new GeneralComparator("at.tugraz.genome.lda.vos.DoubleStringVO", "getValue", "java.lang.Double"));
      
      int row = 1;
      for (int i=(nameAreaVO.size()-1); i!=-1; i--){
        DoubleStringVO vo = nameAreaVO.get(i);
        label = new JLabel(vo.getKey());
        centerPanel.add(label,new GridBagConstraints(0, row, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(1, 3, 1, 3), 0, 0));
        double value = vo.getValue();
        float relArea = area*(float)value;
        if (nameAreaVO.size()==1){
          value = 1d;
          relArea = area;
        }        
        label = new JLabel(Calculator.FormatNumberToString(value*100d, 2)+"%");
        centerPanel.add(label,new GridBagConstraints(1, row, 1, 1, 0.0, 0.0
            ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(1, 3, 1, 3), 0, 0));
        label = new JLabel(String.valueOf(relArea));
        centerPanel.add(label,new GridBagConstraints(2, row, 1, 1, 0.0, 0.0
            ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(1, 3, 1, 3), 0, 0));        
        row++;
      }
      label = new JLabel("Please click accept when you want to store this recommendation");
      centerPanel.add(label,new GridBagConstraints(0, row, 3, 1, 0.0, 0.0
          ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(1, 1, 1, 1), 0, 0));
    }
    
    this.add(centerPanel,BorderLayout.CENTER);
  }
  
  /**
   * initializes the panel containing the buttons "Accept" and "Cancel"
   */
  private void initButtonPanel(){
    JPanel buttonPanel = new JPanel();
    if (acceptAllowed_){
      JButton okButton = new JButton("Accept");
      okButton.setActionCommand("AcceptMSnRecalculation");
      okButton.setToolTipText(TooltipTexts.ACCEPT_GENERAL);
      buttonPanel.add(okButton);
      okButton.addActionListener(this);
      okButton.addActionListener(parent_);
    }
    JButton cancelButton = new JButton("Cancel");
    cancelButton.setActionCommand("DeclineMSnRecalculation");
    cancelButton.setToolTipText(TooltipTexts.CANCEL_GENERAL);
    buttonPanel.add(cancelButton);
    cancelButton.addActionListener(this);
    cancelButton.addActionListener(parent_);
    setDefaultCloseOperation(DO_NOTHING_ON_CLOSE);
    this.add(buttonPanel,BorderLayout.SOUTH);
    
  }

  public void actionPerformed(ActionEvent e)
  {
    if (e.getActionCommand().equalsIgnoreCase("AcceptMSnRecalculation") || 
        e.getActionCommand().equalsIgnoreCase("DeclineMSnRecalculation")){
      setVisible(false); 
      dispose();
    }
  }

  /**
   * 
   * @return the LipidParameterSet containing the new MSn assignments
   */
  public LipidParameterSet getResult()
  {
    return result_;
  }

}
