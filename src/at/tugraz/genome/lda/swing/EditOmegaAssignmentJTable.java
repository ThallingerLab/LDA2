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

package at.tugraz.genome.lda.swing;

import java.awt.BorderLayout;
import java.awt.Container;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Vector;

import javax.swing.Box;
import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.JTextField;
import javax.swing.ListSelectionModel;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.table.DefaultTableModel;

import org.jogamp.java3d.utils.applet.MainFrame;

import at.tugraz.genome.lda.LipidomicsConstants;
import at.tugraz.genome.lda.TooltipTexts;
import at.tugraz.genome.lda.WarningMessage;
import at.tugraz.genome.lda.exception.NumberOutOfRangeException;
import at.tugraz.genome.lda.msn.vos.FattyAcidVO;
import at.tugraz.genome.lda.quantification.LipidParameterSet;
import at.tugraz.genome.lda.utils.StaticUtils;
import at.tugraz.genome.lda.verifier.DoubleVerifier;
import at.tugraz.genome.lda.verifier.IntegerVerifier;
import at.tugraz.genome.lda.vos.DoubleBondPositionVO;

/**
 * This class displays a dialog for editing omega double bond position data
 * 
 * @author Leonida M. Lamp
 *
 */
public class EditOmegaAssignmentJTable extends JFrame implements ActionListener
{

  private static final long serialVersionUID = 1L;
  
  /** the parent action listener */
  private ActionListener parent_;
  /** the title String */
  private final static String TITLE = "Edit \u03C9 - double bond assignment";
  /** button and respective command Strings */
  private final static String BUTTON_CHANGE_SELECTED_SPECIES = "Change Selected";
  private final static String COMMAND_CHANGE_SELECTED_SPECIES = "changeSelectedSpecies";
  private final static String BUTTON_ADD_NEW_SPECIES = "Add New";
  private final static String COMMAND_ADD_NEW_SPECIES = "addNewSpecies";
  private final static String BUTTON_REMOVE_SELECTED_SPECIES = "Remove Selected";
  private final static String COMMAND_REMOVE_SELECTED_SPECIES = "removeSelectedSpecies";
  private final static String BUTTON_ASSIGN_SELECTED_SPECIES = "Assign Selected";
  private final static String COMMAND_ASSIGN_SELECTED_SPECIES = "assignSelectedSpecies";
  private final static String BUTTON_SAVE_CHANGES = "Save Changes";
  private final static String COMMAND_SAVE_CHANGES = "saveOmegaAssignmentChanges";
  private final static String BUTTON_CANCEL = "Cancel";
  private final static String COMMAND_CANCEL = "cancel";
  /** Strings to indicate the accuracy of retention time matches */
  private final String[] accuracy_ = { "low", "medium", "high" };
  /** Textfields for user input */
  private JTextField positionChainFirst_, positionChainSecond_, positionChainThird_, expectedRT_;
  /** local lipidParameterSet for use within this class and lipidParameterSet for saving the changes */
  private LipidParameterSet localLipidParameterSet_, lipidParameterSet_;
  /** Molecular species String without double bond positions of the species being edited */
  private String molecularSpecies_;
  /** Vector of FattyAcidVOs the molecular species being edited is composed of */
  private Vector<FattyAcidVO> chainCombination_; 
  /** local copy of the Vector of DoubleBondPositionVOs storing the omega double bond position information  */
  private Vector<DoubleBondPositionVO> localDoubleBondPositionVOs_;
  /** JTable to display options for omega double bond data */
  private JTable table_;
  /** DefaultTableModel for table_ */
  private DefaultTableModel model_;
  /** static boolean to be false if an error was reported upon trying to save changes */
  private static boolean saveLipidParameterSet_;
  
  
  /**
   * Constructor for a dialog window to edit omega double bond position data
   * @param param LipidParameterSet for which omega double bond position data should be edited
   * @param molecularSpecies Molecular species String without double bond positions of the species being edited
   * @param chainCombination Vector of FattyAcidVOs the species being edited is composed of
   * @param parent the parent ActionListener
   * @param frame the MainFrame for positioning
   */
  public EditOmegaAssignmentJTable(LipidParameterSet param, String molecularSpecies, Vector<FattyAcidVO> chainCombination, ActionListener parent, MainFrame frame) 
  {
    super(TITLE);
    
    this.lipidParameterSet_ = param;
    this.localLipidParameterSet_ = new LipidParameterSet(param); 
    this.molecularSpecies_ = molecularSpecies;
    this.localDoubleBondPositionVOs_ = StaticUtils.getDoubleBondAssignmentsOfMolecularSpecies(
        localLipidParameterSet_.getOmegaInformation(), molecularSpecies_);
    this.chainCombination_ = chainCombination;
    this.parent_ = parent;
    saveLipidParameterSet_ = true;
    
    //general settings
    setLocationRelativeTo(frame);
    setDefaultCloseOperation(DISPOSE_ON_CLOSE);
    setSize(850, 300);
    setVisible(true);
    
    //create panel located north
    String[] individualChains = getIndividualChains(molecularSpecies_);
    JPanel panelNorth = createPanelNorth(individualChains);
    
    //create table located center
    String[] headers = { "Molecular Species", "RT Deviation /s", "Expected RT /min", "Accuracy", "Assigned" };
    
    model_ = new DefaultTableModel() 
    {
      private static final long serialVersionUID = 1L;
      
      @Override
      public boolean isCellEditable(int row, int column) 
      {
        return false;
      }
      
      @Override
      public Class<?> getColumnClass(int column) 
      {
        switch (column) {
          case 4:
            return Boolean.class;
          default:
            return String.class;
        }
      }
    };
    
    for (String header : headers) 
    {
      model_.addColumn(header);
    }
    addRows();
    
    table_ = new JTable(model_);
    table_.setAutoCreateRowSorter(false);
    
    table_.getSelectionModel().addListSelectionListener(new ListSelectionListener()
    {
      public void valueChanged(ListSelectionEvent event) 
      {
        ListSelectionModel lsm = (ListSelectionModel)event.getSource();
        if (lsm.isSelectionEmpty()) 
        {
          //do nothing
        } 
        else 
        {
          int selectedRow = table_.getSelectedRow();
          DoubleBondPositionVO labeledChainCombinationVO = localDoubleBondPositionVOs_.get(selectedRow);
          Vector<FattyAcidVO> fattyAcidVOs = labeledChainCombinationVO.getChainCombination();
          int count = 0;
          expectedRT_.setText(getStringWithTwoDecimals(labeledChainCombinationVO.getExpectedRetentionTime()));
          for (FattyAcidVO fattyAcidVO : fattyAcidVOs) 
          {
            if (chainCombination_.get(count).getDoubleBonds() > 0) 
            {
              String value = String.valueOf(fattyAcidVO.getOmegaPosition());
              if (value.equals("-1")) 
              {
                value=null;
              }
              switch(count) 
              {
                case 0:
                  positionChainFirst_.setText(value);
                  break;
                case 1:
                  positionChainSecond_.setText(value);
                  break;
                case 2:
                  positionChainThird_.setText(value);
                  break;
              }
            }
            count++;
          }
        }
      }
    });
    
    //create panel located south
    JPanel panelSouth = createPanelSouth();
    
    //create container
    Container container = getContentPane();
    container.add(panelNorth, BorderLayout.NORTH);
    container.add(new JScrollPane(table_), BorderLayout.CENTER);
    container.add(panelSouth, BorderLayout.SOUTH);
  } 
  
  /**
   * @return command to save changes
   */
  public static String getSaveChangesCommand() 
  {
    return COMMAND_SAVE_CHANGES;
  }
  
  /**
   * Converts a float to a String with two decimals
   * @param floatValue the float Object
   * @return String with two decimals
   */
  private String getStringWithTwoDecimals(float floatValue) 
  {
    return String.format("%.2f", floatValue);
  }
    
  /**
   * @param molecularSpecies the molecular species
   * @return Strings of individual chains a molecular species is composed of
   */
  private String[] getIndividualChains(String molecularSpecies) 
  {
    String[] individualChains = new String[] {molecularSpecies};
    if (molecularSpecies.contains(LipidomicsConstants.CHAIN_SEPARATOR_NO_POS)) 
    {
      individualChains = molecularSpecies.split(LipidomicsConstants.CHAIN_SEPARATOR_NO_POS);
    } 
    else if (molecularSpecies.contains(LipidomicsConstants.CHAIN_SEPARATOR_KNOWN_POS)) 
    {
      individualChains = molecularSpecies.split(LipidomicsConstants.CHAIN_SEPARATOR_KNOWN_POS);
    }
    return individualChains;
  }
    
  /**
   * Creates a new JButton Object
   * @param buttonText the text to be displayed on the button
   * @param commandText the action command text
   * @param toolTipText the tool tip text
   * @return JButton Object
   */
  private JButton createJButton(String buttonText, String commandText, String toolTipText) 
  {
    JButton jButton = new JButton(buttonText);
    jButton.setActionCommand(commandText);
    jButton.setToolTipText(toolTipText);
    jButton.addActionListener(this);
    return jButton;
  }
    
  /**
   * Invoked when an action occurs
   */
  public void actionPerformed(ActionEvent e)
  {
    int selectedRow = table_.getSelectedRow();
    int[] omegaPositions = new int[]{-1,-1,-1};
    float expectedRT = getMeasuredRT();
    
    switch (e.getActionCommand()) 
    {
      case COMMAND_CHANGE_SELECTED_SPECIES:
        if (selectedRow > -1) 
        {
          DoubleBondPositionVO selected = localDoubleBondPositionVOs_.get(selectedRow);
          float currentRT = selected.getExpectedRetentionTime();
          try 
          {
            expectedRT = Float.parseFloat(expectedRT_.getText());
            selected.setExpectedRetentionTime(expectedRT);
            determineAccuracy(selected);
          } 
          catch (NumberOutOfRangeException ex) 
          {
            int dialogResult = getConfirmDialogResult(expectedRT);
            if (dialogResult == JOptionPane.NO_OPTION || 
                dialogResult == JOptionPane.CLOSED_OPTION) {
              selected.setExpectedRetentionTime(currentRT);
              break;
            }
          } 
          catch (Exception ignore) 
          {
          }
          Vector<FattyAcidVO> fattyAcidVOs = selected.getChainCombination();
          int count = 0;
          for (FattyAcidVO fattyAcidVO : fattyAcidVOs) 
          {
            if (chainCombination_.get(count).getDoubleBonds() > 0) 
            {
              switch(count) 
              {
                case 0:
                  try 
                  {
                    omegaPositions[count] = Integer.parseInt(positionChainFirst_.getText()); 
                    fattyAcidVO.setOmegaPosition(omegaPositions[count]);
                  } 
                  catch (Exception ignore) {}
                  break;
                case 1:
                  try 
                  {
                    omegaPositions[count] = Integer.parseInt(positionChainSecond_.getText()); 
                    fattyAcidVO.setOmegaPosition(omegaPositions[count]);
                  } 
                  catch (Exception ignore) {}
                  break;
                case 2:
                  try 
                  {
                    omegaPositions[count] = Integer.parseInt(positionChainThird_.getText()); 
                    fattyAcidVO.setOmegaPosition(omegaPositions[count]);
                  } 
                  catch (Exception ignore) {}
                  break;
              }
            }
            count++;
          }
          
          model_.setValueAt(selected.getDoubleBondPositionsHumanReadable(), selectedRow, 0);
          model_.setValueAt(getStringWithTwoDecimals(getDeltaRT(expectedRT)), selectedRow, 1);
          model_.setValueAt(expectedRT, selectedRow, 2);
          model_.setValueAt(accuracy_[selected.getAccuracy()], selectedRow, 3);
          model_.setValueAt(selected.getIsAssigned(), selectedRow, 4);
        } 
        else 
        {
          new WarningMessage(new JFrame(), "Error", "No element was selected.");
        }
        break;
      
      case COMMAND_ADD_NEW_SPECIES:
        try {omegaPositions[0] = Integer.parseInt(positionChainFirst_.getText()); } catch (Exception ignore) {}
        try {omegaPositions[1] = Integer.parseInt(positionChainSecond_.getText()); } catch (Exception ignore) {}
        try {omegaPositions[2] = Integer.parseInt(positionChainThird_.getText()); } catch (Exception ignore) {}
        try {expectedRT = Float.parseFloat(expectedRT_.getText());} catch (Exception ignore) {}
        boolean anyOmegaPositionsSet = false;
        for (int omegaPosition : omegaPositions) 
        {
          if (omegaPosition > 0) 
          {
            anyOmegaPositionsSet = true;
            break;
          }
        }
        if (!anyOmegaPositionsSet) 
        {
          new WarningMessage(new JFrame(), "Error", "Select at least one valid \u03C9 - double bond position!");
        } 
        else 
        {
          DoubleBondPositionVO newDoubleBondPositionVO = initiateNewDoubleBondPositionVO(omegaPositions, expectedRT);
          try 
          {
            determineAccuracy(newDoubleBondPositionVO);
          } 
          catch (NumberOutOfRangeException ex) 
          {
            int dialogResult = getConfirmDialogResult(expectedRT);
            if (dialogResult == JOptionPane.NO_OPTION ||
                dialogResult == JOptionPane.CLOSED_OPTION) {
              break;
            }
          }
          localDoubleBondPositionVOs_.add(newDoubleBondPositionVO);
          localLipidParameterSet_.addOmegaInformation(newDoubleBondPositionVO);
          Object[] newRow = createRow(newDoubleBondPositionVO);
          model_.addRow(newRow);
          table_.setModel(model_);
        }
        break;
      
      case COMMAND_ASSIGN_SELECTED_SPECIES:
        if (selectedRow > -1) 
        {
          for (int i=0; i < table_.getRowCount(); i++) 
          {
            model_.setValueAt(false, i, 4);
          }
          localDoubleBondPositionVOs_.get(selectedRow).setIsAssigned(true);
          model_.setValueAt(true, selectedRow, 4);
        } 
        else 
        {
          new WarningMessage(new JFrame(), "Error", "No element was selected.");
        }
        break;
        
      case COMMAND_REMOVE_SELECTED_SPECIES:
        int[] selectedRows = table_.getSelectedRows();
        table_.clearSelection();
        for(int i=0;i<selectedRows.length;i++)
        {
          DoubleBondPositionVO removed = localDoubleBondPositionVOs_.remove(selectedRow);
          localLipidParameterSet_.getOmegaInformation().remove(removed);
          model_.removeRow(selectedRows[i]-i);
        }
        break;
        
      case COMMAND_SAVE_CHANGES:
        Vector<DoubleBondPositionVO> temp = lipidParameterSet_.getOmegaInformation();
        lipidParameterSet_.setOmegaInformation(localLipidParameterSet_.getOmegaInformation());
        parent_.actionPerformed(e);
        if (!saveLipidParameterSet_) 
        {
          lipidParameterSet_.setOmegaInformation(temp);
        } 
        else 
        {
          dispose(); 
        }
        break;
        
      case COMMAND_CANCEL:
        setVisible(false); 
        dispose(); 
        break;
    }
  }
    
    /**
     * Opens a JOptionPane dialog confirming whether the user would like to add an entry with an expected retention time outside the peak range
     * @param expectedRT the expected retention time
     * @return the result of the dialog
     */
    private int getConfirmDialogResult(float expectedRT) 
    {
      int dialogResult = JOptionPane.showConfirmDialog(null, 
          String.format(
              "<html>The expected retention time (%s min) is not inside the range of this peak!<br>"
              + "(Peak range: %s - %s min)<br>"
              + "Would you like to keep these changes anyway?</html>", 
              expectedRT, 
              getStringWithTwoDecimals(localLipidParameterSet_.getIsotopicProbes().get(0).get(0).LowerValley/60f), 
              getStringWithTwoDecimals(localLipidParameterSet_.getIsotopicProbes().get(0).get(0).UpperValley/60f)),
          "Warning",JOptionPane.YES_NO_OPTION);
      return dialogResult;
    }
    
    /**
     * Setter for the static boolean saveLipidParameterSet_, which should be set to false when an error occurred upon trying to save changes
     * @param saveLipidParameterSet the boolean
     */
    public static void setSaveLipidParameterSet(boolean saveLipidParameterSet) 
    {
      saveLipidParameterSet_ = saveLipidParameterSet;
    }
    
    /**
     * Initiates a new DoubleBondPositionVO Object based on the Vector of FattyAcidVOs the molecular species being edited is composed of
     * @param omegaPositions double bond positions (in omega annotation) of the individual FattyAcidVOs
     * @param expectedRT the expected retention time 
     * @return newly created DoubleBondPositionVO object
     */
    private DoubleBondPositionVO initiateNewDoubleBondPositionVO(int[] omegaPositions, float expectedRT)
    {
      Vector<FattyAcidVO> chainCombination = new Vector<FattyAcidVO>();
      int count = 0;
      for (FattyAcidVO fattyAcid : chainCombination_) {
        chainCombination.add(new FattyAcidVO(fattyAcid));
        chainCombination.get(count).setOmegaPosition(omegaPositions[count]);
        count++;
      }
      
      DoubleBondPositionVO newDoubleBondPositionVO = new DoubleBondPositionVO(
          chainCombination, expectedRT, DoubleBondPositionVO.ACCURACY_LOW, molecularSpecies_, false);
      return newDoubleBondPositionVO;
    }
    
    /**
     * Determines the accuracy of a double bond position based on its expected retention time
     * @param doubleBondPositionVO the double bond position value object
     * @throws NumberOutOfRangeException when the given expected retention time is outside the peak range
     */
    private void determineAccuracy(DoubleBondPositionVO doubleBondPositionVO) throws NumberOutOfRangeException 
    {
      float expectedRT = doubleBondPositionVO.getExpectedRetentionTime();
      Range[] peakRanges = StaticUtils.determinePeakRanges(localLipidParameterSet_);
      Range peakLimits = peakRanges[0];
      Range mediumAccuracy = peakRanges[1];
      Range highAccuracy = peakRanges[2];
      
      doubleBondPositionVO.setAccuracy(DoubleBondPositionVO.ACCURACY_LOW);
      if (peakLimits.insideRange(expectedRT)) 
      {
        if (mediumAccuracy.insideRange(expectedRT)) doubleBondPositionVO.setAccuracy(DoubleBondPositionVO.ACCURACY_MEDIUM);
        if (highAccuracy.insideRange(expectedRT)) doubleBondPositionVO.setAccuracy(DoubleBondPositionVO.ACCURACY_HIGH);
      } 
      else 
      {
        throw new NumberOutOfRangeException();
      }
    }
    
    /**
     * @return the measured retention time of the lipidParameterSet being edited as float object
     */
    private float getMeasuredRT() 
    {
      return Float.parseFloat(localLipidParameterSet_.getRt());
    }
    
    /**
     * Calculates the time difference of the measured and expected retention times in seconds
     * @param expectedRT the expected retention time
     * @return the retention time difference as float object
     */
    private float getDeltaRT(float expectedRT) 
    {
      return Math.abs(getMeasuredRT() - expectedRT) *60;
    }
    
    /**
     * Creates the JPanel to be located north
     * @param individualChains individual chains the molecular species is composed of
     * @return JPanel object
     */
    private JPanel createPanelNorth(String[] individualChains) 
    {
      JPanel panelNorth = new JPanel(new BorderLayout());
      JPanel addNewSpeciesPanel = new JPanel();
      JPanel measuredRTPanel = new JPanel();
      
      JLabel label = new JLabel(String.format("Measured RT: %s min", getStringWithTwoDecimals(getMeasuredRT())));  
      label.setToolTipText(TooltipTexts.DISPLAY_MEASURED_RT);
      measuredRTPanel.add(label);
      
      JButton changeSelectedSpeciesButton = createJButton(BUTTON_CHANGE_SELECTED_SPECIES, COMMAND_CHANGE_SELECTED_SPECIES, TooltipTexts.DISPLAY_CHANGE_SELECTED_DB_SPECIES);
      JButton addNewSpeciesButton = createJButton(BUTTON_ADD_NEW_SPECIES, COMMAND_ADD_NEW_SPECIES, TooltipTexts.DISPLAY_ADD_NEW_DB_SPECIES);
      
      int count = 0;
      for (String chain : individualChains) 
      {
        if (chainCombination_.get(count).getDoubleBonds() > 0) 
        {
          label = new JLabel(String.format("%s n-", chain));
          label.setToolTipText(TooltipTexts.DISPLAY_CHAIN_POSITION);
          addNewSpeciesPanel.add(label);
          switch(count) 
          {
            case 0:
              positionChainFirst_ = new JTextField(2);
              positionChainFirst_.setInputVerifier(new IntegerVerifier());
              positionChainFirst_.setToolTipText(TooltipTexts.DISPLAY_CHAIN_POSITION);
              addNewSpeciesPanel.add(positionChainFirst_);
              break;
            case 1:
              positionChainSecond_ = new JTextField(2);
              positionChainSecond_.setInputVerifier(new IntegerVerifier());
              positionChainSecond_.setToolTipText(TooltipTexts.DISPLAY_CHAIN_POSITION);
              addNewSpeciesPanel.add(positionChainSecond_);
              break;
            case 2:
              positionChainThird_ = new JTextField(2);
              positionChainThird_.setInputVerifier(new IntegerVerifier());
              positionChainThird_.setToolTipText(TooltipTexts.DISPLAY_CHAIN_POSITION);
              addNewSpeciesPanel.add(positionChainThird_);
              break;
          }
        }
        count++;
      }
      addNewSpeciesPanel.add(Box.createHorizontalStrut(5));
      label = new JLabel("Expected RT /min: "); 
      label.setToolTipText(TooltipTexts.DISPLAY_EXPECTED_RT);
      addNewSpeciesPanel.add(label);
      expectedRT_ = new JTextField(4);
      expectedRT_.setInputVerifier(new DoubleVerifier(true));
      expectedRT_.setToolTipText(TooltipTexts.DISPLAY_EXPECTED_RT);
      addNewSpeciesPanel.add(expectedRT_);
      addNewSpeciesPanel.add(Box.createHorizontalStrut(5));
      addNewSpeciesPanel.add(changeSelectedSpeciesButton);
      addNewSpeciesPanel.add(addNewSpeciesButton);
      
      panelNorth.add(measuredRTPanel, BorderLayout.EAST);
      panelNorth.add(addNewSpeciesPanel, BorderLayout.WEST);
      
      return panelNorth;
   }
   
   /**
    * Creates JPanel located south
    * @return JPanel object
    */
   private JPanel createPanelSouth() 
   {
     JPanel panelSouth = new JPanel(new BorderLayout());
     JPanel editPanel = new JPanel();
     JPanel cancelPanel = new JPanel();
     JButton assignSelectedSpeciesButton = createJButton(BUTTON_ASSIGN_SELECTED_SPECIES, COMMAND_ASSIGN_SELECTED_SPECIES, TooltipTexts.DISPLAY_ASSIGN_SELECTED_DB_SPECIES);
     JButton removeSelectedSpeciesButton = createJButton(BUTTON_REMOVE_SELECTED_SPECIES, COMMAND_REMOVE_SELECTED_SPECIES, TooltipTexts.DISPLAY_REMOVE_SELECTED_DB_SPECIES);
     JButton saveButton = createJButton(BUTTON_SAVE_CHANGES, COMMAND_SAVE_CHANGES, TooltipTexts.DISPLAY_SAVE_CHANGES);
     JButton cancelButton =  createJButton(BUTTON_CANCEL, COMMAND_CANCEL, TooltipTexts.CANCEL_GENERAL);
     editPanel.add(assignSelectedSpeciesButton);
     editPanel.add(removeSelectedSpeciesButton);
     editPanel.add(saveButton);
     cancelPanel.add(cancelButton);
     panelSouth.add(editPanel, BorderLayout.WEST);
     panelSouth.add(cancelPanel, BorderLayout.EAST);
     
     return panelSouth;
   }
   
   /**
    * Adds rows to model_
    */
   private void addRows() 
   {
     for (DoubleBondPositionVO doubleBondPositionVO : localDoubleBondPositionVOs_) 
     {
       Object[] row = createRow(doubleBondPositionVO);
       model_.addRow(row);
     }
   }
   
   /**
    * Creates rows from DoubleBondPositionVO data
    * @param doubleBondPositionVO the double bond position value object
    * @return Array of row Objects
    */
   private Object[] createRow(DoubleBondPositionVO doubleBondPositionVO) 
   {
     Object[] row = { doubleBondPositionVO.getDoubleBondPositionsHumanReadable(), 
         getStringWithTwoDecimals(getDeltaRT(doubleBondPositionVO.getExpectedRetentionTime())),
         getStringWithTwoDecimals(doubleBondPositionVO.getExpectedRetentionTime()),
         accuracy_[doubleBondPositionVO.getAccuracy()], 
         doubleBondPositionVO.getIsAssigned()
     };
     return row;
   }
    

  }
