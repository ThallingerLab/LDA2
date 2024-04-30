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
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Image;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.Vector;

import javax.imageio.ImageIO;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSplitPane;
import javax.swing.JTabbedPane;
import javax.swing.JTextField;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import at.tugraz.genome.lda.interfaces.rdi.GeneralSettingsPanelListener;
import at.tugraz.genome.lda.listeners.AddChainEquationDocumentListener;
import at.tugraz.genome.lda.listeners.AddChainFragmentDocumentListener;
import at.tugraz.genome.lda.listeners.AddHeadEquationDocumentListener;
import at.tugraz.genome.lda.listeners.AddHeadFragmentDocumentListener;
import at.tugraz.genome.lda.listeners.AddPositionEquationDocumentListener;
import at.tugraz.genome.lda.listeners.ChangeEquationDocumentListener;
import at.tugraz.genome.lda.listeners.ChangeFragmentDocumentListener;
import at.tugraz.genome.lda.listeners.ChangeFragmentFocusListener;
import at.tugraz.genome.lda.msn.parser.FragRuleParser;
import at.tugraz.genome.lda.LipidomicsConstants;
import at.tugraz.genome.lda.Settings;
import at.tugraz.genome.lda.WarningMessage;
import at.tugraz.genome.lda.msn.FragmentCalculator;
import at.tugraz.genome.lda.msn.LipidomicsMSnSet;
import at.tugraz.genome.lda.msn.MSnAnalyzer;
import at.tugraz.genome.lda.exception.ChemicalFormulaException;
import at.tugraz.genome.lda.exception.HydroxylationEncodingException;
import at.tugraz.genome.lda.exception.LipidCombinameEncodingException;
import at.tugraz.genome.lda.exception.NoRuleException;
import at.tugraz.genome.lda.exception.RulesException;
import at.tugraz.genome.lda.msn.RulesContainer;
import at.tugraz.genome.lda.msn.vos.FragmentRuleVO;
import at.tugraz.genome.lda.msn.vos.IntensityChainVO;
import at.tugraz.genome.lda.msn.vos.IntensityRuleVO;
import at.tugraz.genome.lda.msn.vos.MSnDebugVO;
import at.tugraz.genome.lda.quantification.LipidParameterSet;
import at.tugraz.genome.lda.quantification.LipidomicsAnalyzer;
import at.tugraz.genome.lda.swing.rdi.FragmentFormulaTextField;
import at.tugraz.genome.lda.swing.rdi.FragmentNameTextField;
import at.tugraz.genome.lda.swing.rdi.GeneralSettingsPanel;
import at.tugraz.genome.lda.swing.rdi.IntegerRangeDocument;
import at.tugraz.genome.lda.utils.StaticUtils;
import at.tugraz.genome.lda.verifier.EquationVerifier;
import at.tugraz.genome.lda.verifier.FragmentVerifier;
import at.tugraz.genome.lda.vos.rdi.GeneralSettingsVO;
import at.tugraz.genome.maspectras.parser.exceptions.SpectrummillParserException;
import at.tugraz.genome.maspectras.parser.spectrummill.ElementConfigParser;
import at.tugraz.genome.maspectras.quantification.CgAreaStatus;
import at.tugraz.genome.maspectras.quantification.CgException;
import at.tugraz.genome.maspectras.quantification.CgProbe;


/**
 * Class builds the interface for new rules
 * @author Andreas Ziegl
 * @author Juergen Hartler
 *
 */
public class RuleDefinitionInterface extends JSplitPane implements GeneralSettingsPanelListener
{  
  private static final long serialVersionUID = 6626019779002700255L;  
   
  /** Split Pane at the top */
  public JSplitPane topSplitPane_;
  
  /** Layout for the hole interface */
  private GridBagLayout gridBayLayout_;
  
  /** Upper Section always visible */
  private JScrollPane permanentHead_ = new JScrollPane(); 
  
  /** The middle section who depend form the class */
  private JTabbedPane middleTabSection_;
  
  /** Bottom Section */
  private JPanel bottomSectionButtons_;
  
  /** Lipid Parameters */
  private LipidParameterSet data_;
  
  private int highestMSLevel_;
  
  /** Name of the lipid class */
  private  String lipidClassName_;  
  
  /** Name of the lipid */
  private  String lipidName_;
  
  /**Cache if there a problems with the spectra */
//  private  String lipidNameBefore_;
  
  /** Name of the adduct of the lipid */
  private  String lipidAdduct_;    
  
  /** Precursor value of the lipid */
  private  double lipidPrecursorValue_;
  
  /** Chemical formula of the lipid */
  private String lipidFormula_;
  
  private GeneralSettingsVO generalSettingsVO_;
  
  /**Cache if there a problems with the spectra */
//  private  double lipidPrecursorValueBefore_;  
  
  /** Name of the class to identify the rules */
  private String ruleClassIdentifier_;
  
  /** Head fragment rules set */
  private java.util.Hashtable<String,FragmentRuleVO> headFragmentRules_;  
  
  /** Head intensity rules set */
  private Vector<IntensityRuleVO> headIntensityRules_;
  
  /** Chain fragment rules set */
  private java.util.Hashtable<String,FragmentRuleVO> chainFragmentRules_;
  
  /** Chain intensity rules set */
  private Vector<IntensityRuleVO> chainIntensityRules_;
  
  /** Position intensity rules set */
  private Vector<IntensityRuleVO> positionIntensityRules_;
  
  
  /** Sorted head fragment rules */
  private Vector<FragmentRuleVO> headFragments_;
  
  /** Sorted chain fragment rules - save rules accessed to this vector */
  private Vector<FragmentRuleVO> sortedChainFragments_;
  
  /** the element parser */
  private ElementConfigParser elementParser_ = Settings.getElementParser();  
  
  /** Possibility List for combo mandatory fields */
  private String possibibilityListPostProcessing_[] = {"true", "false"};
      
  /** The analyzer from LipidomicsAnalyzer as a parameter of this class*/
  LipidomicsAnalyzer analyzer_;
  
  /** Field for a new fragment name */
  private JTextField fragRuleNameField_;
  
  /** Field for a new fragment formula */
  private JTextField fragRuleFormulaField_; 
  
  /** Field for a new fragment charge */
  private JTextField fragRuleChargeField_; 
  
  /** Field for a new fragment ms level */
  private JTextField fragRuleMsLevelField_;
  
  /** Field for a new chainfragment name */
  private JTextField chainRuleNameField_;
  
  /** Field for a new chainfragment formula */
  private JTextField chainRuleFormulaField_;
  
  /** Field for a new chainfragment charge */
  private JTextField chainRuleChargeField_;
  
  /** Field for a new chainfragment ms level */
  private JTextField chainRuleMsLevelField_;
  
  /** Field for a new fragment mandatory */
  private JComboBox<Object> fragRuleMandatoryCombo_;
  
  /** Field for a new chainfragment mandatory */
  private JComboBox<Object> chainRuleMandatoryCombo_;
  
  /** Sign to check the head fragment */
  private JLabel fragRuleOK_;
  
  /** Sign to check the chain fragment */
  private JLabel chainRuleOK_;
  
  /** Possible Sign graphic green */
  private java.net.URL resource1_;
  
  /** Possible Sign graphic red */
  private java.net.URL redCrossPicture_;
  
  /** Tab for the head rules */
  private JPanel headRuleTab_;
  
  /** Tab for the chain rules */
  private JPanel chainTab_;
  
  /** Tab for the position rules */
  private JPanel positionTab_;
  
  /** Field to enter a new head fragment equation */
  private JTextField equRuleEquationField_;
  
  /** Field to enter a new head fragment mandatory */
  private JComboBox<Object> equRuleMandatoryCombo_;
  
  /** Sign to check the head fragment */
  private JLabel equRuleOK_;
  
  /** Sign to check the position equation */  
  private JLabel positionRuleOK_;
  
  /** Field to enter a new position equation */
  private JTextField positionEquationField_;
  
  /** Field to choose a new mandatory for a position equation */
  private JComboBox<Object> positionRuleMandatoryCombo_;
  
  /** Sign graphic space if the new rule is ok */
  private JLabel equChainRuleOK_;
  
  /** Field to enter a new chain equation */
  private JTextField equChainRuleEquationField_;
  
  /** Field to choose a new mandatory for a chain equation */
  private JComboBox<Object> equChainRuleMandatoryCombo_;
  
  /** The main object of the software */
  private SpectrumUpdateListener spectrumUpdater_;
  
  /** The show details box */
  private JDialog showHeadDetails_;
  
  /** The global msn Analyzer */
  MSnAnalyzer msnAnalyzer_;
  
  /** Fragments are only checked when errorMessageNeeded is true */
  public boolean errorMessageNeeded_ = false;
  
  /** Verifier for the formulas of head fragments */
  public FragmentVerifier[] headFragmentFormulaVerifier_;
  
  /** Verifier for the charges of head fragments */
  public FragmentVerifier[] headFragmentChargeVerifier_;
  
  /** Verifier for the mslevels of head fragments */
  public FragmentVerifier[] headFragmentMsLevelVerifier_;  

  /** the name of each head fragment */
  private JTextField headFragRuleName_;
  
  /** The head rule formula textfields */
  public static JTextField[] headRuleFormulas_; 
  
  /** The head rule charge textfields */
  public static JTextField[] headRuleCharges_; 
  
  /** The head rule mslevel textfields */
  public static JTextField[] headRuleMSLevels_; 
  
  /** Verifier for the formulas of chain fragments */
  public FragmentVerifier[] chainFragmentFormulaVerifier_;
  
  /** Verifier for the charges of chain fragments */
  public FragmentVerifier[] chainFragmentChargeVerifier_;
  
  /** Verifier for the mslevels of chain fragments */
  public FragmentVerifier[] chainFragmentMsLevelVerifier_;  

  /** the name of each chain fragment */
  private JTextField chainFragRuleName_;
  
  /** The chain rule formula textfields */
  public static JTextField[] chainRuleFormulas_; 
  
  /** The chain rule charge textfields */
  public static JTextField[] chainRuleCharges_; 
  
  /** The chain rule mslevel textfields */
  public static JTextField[] chainRuleMSLevels_;  
  
  /** Public head Fragment combo */
  public JComboBox<Object> fragRuleMandatoryField_;

  /** Head fragment mandatories */	
	private boolean[] headFragmentMandatories_;
	
	/** Chain fragment mandatories */	
	private boolean[] chainFragmentMandatories_;
	
	/** The verifier for head equations */
	public EquationVerifier[] headEquationVerifier_; 
	
	/** The fields for head rules */
  public JTextField[] headRuleEquationFields_;
  
  /** The mandatories of the head equations */
  private boolean[] headEquationMandatories_;
  
  /** The verifiers of the chain equations */
  public EquationVerifier[] chainEquationVerifier_;
  
  /** The equations fields of the chain rules */
  public JTextField[] chainRuleEquationFields_;
  
  /** The mandatories of the chain equations */
  private boolean[] chainEquationMandatories_;
  
  /** The verifiers of the position equations */
  public EquationVerifier[] positionEquationVerifier_;
  
  /** The fields of the position equations */
  public JTextField[] positionRuleEquationFields_;
  
  
  /** The mandatories of the position euqations */
  private boolean[] positionEquationMandatories_;

  /** directory for caching the intermediate rules*/
  public final static String CACHE_DIR = "cache";
    
  /** the currently affected section is head group*/
  private final static int TYPE_HEAD = 1;
  /** the currently affected section is chain*/
  private final static int TYPE_CHAIN = 2;
  /** the currently affected section is position*/
  private final static int TYPE_POSITION = 3;
  
  private final int MIN_MSLEVEL = 2;
  
  /** a VO containing the entered values for the required general settings - works as cache if wrong values are enterd*/
  private GeneralSettingsPanel generalSettings_;
  
  /** true if the tabs for the fragment definition shall be displayed */
  private boolean showFragmentTabs_;

  
/**
 * Creates the hole panel by including the different sections except the bottom section
 */
  public void createPanel()
  {    
    generalTopSection();    
    ruleTabsSection();     
    topSplitPane_ = new JSplitPane(JSplitPane.VERTICAL_SPLIT);    
    topSplitPane_.setDividerSize(0);    
    topSplitPane_.setTopComponent(permanentHead_);    
    topSplitPane_.setBottomComponent(middleTabSection_);    
    topSplitPane_.setResizeWeight(0.4);   
    setTopComponent(topSplitPane_);
//    paintNewSpectra(false);   
  }
  
  /**
   * Creates the bottom section to add it afterward
   * @return bottomSectionButtons_ the section added to the panel
   */
  public JPanel makeFinalButtonSection()
  {
    finalButtonsSection();
    return bottomSectionButtons_;
  }  

  /**
   * Creates the top section with the first general informations
   */
  public void generalTopSection()
  {
    JPanel permanent = new JPanel();
    gridBayLayout_ = new GridBagLayout();
    permanent.setLayout(gridBayLayout_);
    
    JLabel spaceClass1 = new JLabel( "    " );
    gridBayLayout_.setConstraints(spaceClass1, new GridBagConstraints (0, 0, 1, 1, 0, 0, 
    GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0 ));    
    
    JLabel lipidClass = new JLabel( "  Lipid Class  " );
    gridBayLayout_.setConstraints(lipidClass, new GridBagConstraints (0, 1, 1, 1, 0, 0, 
    GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0 ));            
    final JTextField lipidClassField = new JTextField();
    gridBayLayout_.setConstraints(lipidClassField, new GridBagConstraints (1, 1, 1, 1, 1000, 0, 
    GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0 ));
    lipidClassField.setColumns( 15 );  
    lipidClassField.setText(lipidClassName_);  
    lipidClassField.setToolTipText("Enter the lipid class here!");
   
    JLabel adductClass = new JLabel( "  Adduct  " );
    gridBayLayout_.setConstraints(adductClass, new GridBagConstraints (0, 2, 1, 1, 0, 0, 
    GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0 ));      
    final JTextField adductClassField = new JTextField();
    gridBayLayout_.setConstraints(adductClassField, new GridBagConstraints (1, 2, 1, 1, 1000, 0.0, 
    GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0 ));
    adductClassField.setColumns( 15 );     
    adductClassField.setText(lipidAdduct_);  
    adductClassField.setToolTipText("Enter the adduct here!");
    
    JLabel nameClass = new JLabel( "  Name  " );
    gridBayLayout_.setConstraints(nameClass, new GridBagConstraints (0, 3, 1, 1, 0, 0, 
    GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0 ));      
    final JTextField nameClassField = new JTextField();
    gridBayLayout_.setConstraints(nameClassField, new GridBagConstraints (1, 3, 1, 1, 1000, 0.0, 
    GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0 ));
    nameClassField.setColumns( 15 ); 
    nameClassField.setText(lipidName_);
    nameClassField.setEditable( true );
    nameClassField.setToolTipText("Enter the full name with the _RT here!");

    
    JLabel precursorClass = new JLabel( "  Precursor  " );
    gridBayLayout_.setConstraints(precursorClass, new GridBagConstraints (2, 1, 1, 1, 0, 0, 
    GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0 ));      
    final JTextField precursorClassField = new JTextField();
    gridBayLayout_.setConstraints(precursorClassField, new GridBagConstraints (3, 1, 1, 1, 1000, 0.0, 
    GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0 ));
    precursorClassField.setColumns( 15 );  
    precursorClassField.setText(Double.toString(lipidPrecursorValue_));
    precursorClassField.setEditable( true );
    precursorClassField.setToolTipText("Enter the precursor here!");

    
    JLabel formulaLabel = new JLabel( "  Neutral formula  " );
    gridBayLayout_.setConstraints(formulaLabel, new GridBagConstraints (2, 2, 1, 1, 0, 0, 
    GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0 ));      
    final JTextField formulaField = new JTextField();
    gridBayLayout_.setConstraints(formulaField, new GridBagConstraints (3, 2, 1, 1, 1000, 0.0, 
    GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0 ));
    formulaField.setColumns( 15 );  
    formulaField.setText(lipidFormula_);
    formulaField.setEditable( true );
    formulaField.setToolTipText("Enter the precursor here!");

    
    JLabel spaceClass2 = new JLabel( "    " );
    gridBayLayout_.setConstraints(spaceClass2, new GridBagConstraints (0, 4, 1, 1, 0, 0, 
    GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0 ));    
    
    JButton changeSettingsClass = new JButton( "Change Settings" );
    gridBayLayout_.setConstraints(changeSettingsClass, new GridBagConstraints (0, 5, 1, 1, 0, 0, 
    GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0 ));   
    changeSettingsClass.setToolTipText("Changes the settings and refreshes the spectrum");
    changeSettingsClass.addActionListener(new ActionListener() 
    { 
        public void actionPerformed(ActionEvent e)
        {   
          deleteDetailsBoxes();
          boolean error = false;
          ruleClassIdentifier_ = StaticUtils.getRuleName(lipidClassField.getText(), adductClassField.getText()); 
          try
          {       
          getRules();
          }
          catch (NoRuleException e1) 
          {            
            int reply = JOptionPane.showConfirmDialog(topSplitPane_, "This class does not exist. Create new one?", "Create new class" , JOptionPane.YES_NO_OPTION);
            if (reply == JOptionPane.YES_OPTION) 
            { 
              error = false;
              setEverythingNull();
              lipidClassName_ = lipidClassField.getText();
              lipidAdduct_ = adductClassField.getText();
              showFragmentTabs_ = false;
            }           
            else
            {              
              error = true; 
            }
          }
          catch (SpectrummillParserException e1) 
          {            
            JOptionPane.showMessageDialog(topSplitPane_, e1.getMessage()); 
          }
          catch (RulesException e1) 
          {            
            JOptionPane.showMessageDialog(topSplitPane_, e1.getMessage()); 
          }
          catch (IOException e1) 
          {           
            JOptionPane.showMessageDialog(topSplitPane_, e1.getMessage());               
          }
          
          if(error == false)
          {
            lipidClassName_ = lipidClassField.getText();
            lipidAdduct_ = adductClassField.getText();
            try{
              lipidPrecursorValue_ = Double.parseDouble(precursorClassField.getText());  
              lipidFormula_ = formulaField.getText();
              lipidName_ = nameClassField.getText();
            
              data_.setNameString(nameClassField.getText());
            
              Float mz = Float.parseFloat(precursorClassField.getText());
              if (StaticUtils.checkChemicalFormula(formulaField.getText())){
                data_ = new LipidParameterSet(mz, data_.getName(), data_.getDoubleBonds(), adductClassField.getText(),
                  data_.getPreciseRT(), formulaField.getText(), "",1,data_.getOhNumber());
                data_.LowerMzBand = LipidomicsConstants.getCoarseChromMzTolerance(mz);
                data_.UpperMzBand = LipidomicsConstants.getCoarseChromMzTolerance(mz);
            
                refreshMiddleWithoutCurrentGenerals(0);
                paintNewSpectra(true);
              }
            }catch(NumberFormatException nfx){
              new WarningMessage(new JFrame(), "Error", "The entered precursor m/z value is not double format!");
            }
          }          
        }
    }); 
    
    permanent.add(spaceClass1);
    permanent.add(lipidClass);
    permanent.add(lipidClassField);    
    permanent.add(adductClass);
    permanent.add(adductClassField);    
    permanent.add(precursorClass);
    permanent.add(precursorClassField);
    permanent.add(formulaLabel);
    permanent.add(formulaField);
    permanent.add(nameClass);
    permanent.add(nameClassField);    
    permanent.add(spaceClass2);
    permanent.add(changeSettingsClass);      
    permanentHead_ = new JScrollPane(permanent); 
    
  }
  
  /**
   * Creates the four tabs in the middle section to enter new rules
   */
  public void ruleTabsSection()
  {
    middleTabSection_ = new JTabbedPane
    (JTabbedPane.TOP,JTabbedPane.SCROLL_TAB_LAYOUT );  
    createGeneralTab();
    createHeadRuleTab();
    createChainTab();
    createPositionTab();
    middleTabSection_.addChangeListener(new ChangeListener() {
    public void stateChanged(ChangeEvent e) {
        if (showFragmentTabs_ && !generalSettings_.updateGeneralEntries(false)) {
          middleTabSection_.setSelectedIndex(0);         
      }

    }
});
  }
  
   /**
    * Creates the general tab and adds it to the middle tab section 
    */
  public void createGeneralTab()
  {
    generalSettings_ = new GeneralSettingsPanel(this.generalSettingsVO_,this);        
    JScrollPane general = new JScrollPane(generalSettings_);
    middleTabSection_.addTab("General", general);    
  }
   
  
  
  /**
   * Creates the head rule tab and adds it to the middle tab section 
   */
  public void createHeadRuleTab()
  {  		
    headRuleTab_ = new JPanel();    
    headRuleTab_.setLayout(gridBayLayout_);  
    
    resource1_ = RuleDefinitionInterface.class.getResource( "/images/uploaded.gif" );   
    redCrossPicture_ = RuleDefinitionInterface.class.getResource( "/images/error.gif" );    

    JLabel fragRuleName = new JLabel( "Name" );        
    gridBayLayout_.setConstraints(fragRuleName, new GridBagConstraints (1, 0, 1, 1, 0, 0, 
    GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0 ));
   
    JLabel fragRuleFormula = new JLabel( "Formula" );        
    gridBayLayout_.setConstraints(fragRuleFormula, new GridBagConstraints (2, 0, 1, 1, 0, 0, 
    GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0 ));
   
    JLabel fragRuleCharge = new JLabel( "Charge" );        
    gridBayLayout_.setConstraints(fragRuleCharge, new GridBagConstraints (3, 0, 1, 1, 0, 0, 
    GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0 ));
   
    JLabel fragRuleMsLevel = new JLabel( "MS Level" );        
    gridBayLayout_.setConstraints(fragRuleMsLevel, new GridBagConstraints (4, 0, 1, 1, 0, 0, 
    GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0 ));
   
    JLabel fragRuleMandatory = new JLabel( "Mandatory" );        
    gridBayLayout_.setConstraints(fragRuleMandatory, new GridBagConstraints (5, 0, 1, 1, 0, 0, 
    GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0 ));   

      
    
    headFragments_ = getSortedVector(headFragmentRules_, TYPE_HEAD);
    if(headFragments_.size() !=0)
    {
      headFragmentFormulaVerifier_ = new FragmentVerifier [headFragments_.size()];
      headRuleFormulas_ = new JTextField [headFragments_.size()];
      
      headFragmentChargeVerifier_ = new FragmentVerifier [headFragments_.size()];
      headRuleCharges_ = new JTextField [headFragments_.size()];
      
      headFragmentMsLevelVerifier_ = new FragmentVerifier [headFragments_.size()];
      headRuleMSLevels_ = new JTextField [headFragments_.size()];
      
      headFragmentMandatories_ = new boolean [headFragments_.size()];
      
      for(int i=0; i<headFragments_.size(); i++)
      {
    	final int positionInTable = i;  
    	  
        headFragRuleName_ = new JTextField();
        gridBayLayout_.setConstraints(headFragRuleName_, new GridBagConstraints (1, i+1, 1, 1, 0, 0, 
        GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0 ));
        headFragRuleName_.setColumns( 20 );
        headFragRuleName_.setText( headFragments_.elementAt(i).getName() );
        headFragRuleName_.setEditable( false );
       
        FragmentVerifier fragmentFormulaVerifier = new FragmentVerifier(this, FragmentVerifier.TYPE_FORMULA, 1);        
        JTextField headRuleFormula = fragmentFormulaVerifier.getTextField();
        gridBayLayout_.setConstraints(headRuleFormula, new GridBagConstraints (2, i+1, 1, 1, 0, 0, 
        GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0 ));                    
        headRuleFormula.setText( headFragments_.elementAt(i).getFormula() );        
        headRuleFormula.getDocument().addDocumentListener(new ChangeFragmentDocumentListener(this, positionInTable, 0, 1));
        headRuleFormula.addFocusListener(new ChangeFragmentFocusListener(this, positionInTable, 0, 1));

        headFragmentFormulaVerifier_[i] = fragmentFormulaVerifier;
        headRuleFormulas_[i] = headRuleFormula;
        
        FragmentVerifier fragmentChargeVerifier = new FragmentVerifier(this, 2, 1);
        JTextField headRuleCharge = fragmentChargeVerifier.getTextField();
        headRuleCharge.setColumns( 5 );               
        gridBayLayout_.setConstraints(headRuleCharge, new GridBagConstraints (3, i+1, 1, 1, 0, 0, 
        GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0 ));  
        headRuleCharge.setText( Integer.toString(headFragments_.elementAt(i).getCharge()) ); 
        headRuleCharge.setEditable( true );        
        headRuleCharge.getDocument().addDocumentListener(new ChangeFragmentDocumentListener(this, positionInTable, 1, 1));
        headRuleCharge.addFocusListener(new ChangeFragmentFocusListener(this, positionInTable, 1, 1));      
        headFragmentChargeVerifier_[i] = fragmentChargeVerifier;
        headRuleCharges_[i] = headRuleCharge;
        
        
        FragmentVerifier fragmentMsLevelVerifier = new FragmentVerifier(this, 3, 1);
        JTextField headRuleMSLevel = fragmentMsLevelVerifier.getTextField();
        headRuleMSLevel.setDocument(new IntegerRangeDocument(MIN_MSLEVEL,highestMSLevel_));

        headRuleMSLevel.setColumns( 5 );
        headRuleMSLevel.setText( Integer.toString(headFragments_.elementAt(i).getMsLevel()) );
        gridBayLayout_.setConstraints(headRuleMSLevel, new GridBagConstraints (4, i+1, 1, 1, 0, 0, 
        GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0 ));
        headRuleMSLevel.setEditable( true );       
        headRuleMSLevel.getDocument().addDocumentListener(new ChangeFragmentDocumentListener(this, positionInTable, 2, 1));        
        headRuleMSLevel.addFocusListener(new ChangeFragmentFocusListener(this, positionInTable, 2, 1));      
        headFragmentMsLevelVerifier_[i] = fragmentMsLevelVerifier;
        headRuleMSLevels_[i] = headRuleMSLevel;
        
        JComboBox<Object> fragRuleMandatoryField = new JComboBox<Object>(possibibilityListPostProcessing_);          
        if(headFragments_.elementAt(i).isMandatory()==FragmentRuleVO.MANDATORY_TRUE)
        {
          fragRuleMandatoryField.setSelectedIndex(0); 
        }           
        else
        {
          fragRuleMandatoryField.setSelectedIndex(1);
        }          
        gridBayLayout_.setConstraints(fragRuleMandatoryField, new GridBagConstraints (5, i+1, 1, 1, 0, 0, 
        GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0 ));   
        fragRuleMandatoryField.addActionListener(new ActionListener() 
        {
          public void actionPerformed(ActionEvent eventSource) 
          {           	
          	headFragmentMandatories_[positionInTable] = intToBoolean(fragRuleMandatoryField.getSelectedIndex());
          	refreshTextfieldToVerifie(positionInTable, 0, 1);
          	checkFragmentWithErrors(positionInTable, 2, 1);
          }
        });           
              
        JButton deleteButton = new JButton();
        gridBayLayout_.setConstraints(deleteButton, new GridBagConstraints (6, i+1, 1, 1, 0, 0, 
        GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0 ));   
        deleteButton.setToolTipText("Deletes this fragment");
      
        try
        {
          Image img = ImageIO.read(getClass().getResource("/images/error.gif"));
          deleteButton.setIcon(new ImageIcon(img));
        }
        catch (IOException e)
        {      
        }
        final String key = headFragments_.elementAt(i).getName();
        deleteButton.addActionListener(new ActionListener() 
        {
          public void actionPerformed(ActionEvent event) 
          {
            int reply = JOptionPane.showConfirmDialog(topSplitPane_, "Delete \"" + headFragmentRules_.get(key).getName() +
                "\" rule?", "Delete Rule" , JOptionPane.YES_NO_OPTION);
            if (reply == JOptionPane.YES_OPTION) 
            {
              headFragmentRules_.remove(key);
              deleteDependingRules(key,TYPE_HEAD);
              refreshMiddle(1);
              paintNewSpectra(false);                          
            }                
          }
        });
      
        headRuleTab_.add( headFragRuleName_ );
        headRuleTab_.add( headRuleFormulas_[i] );
        headRuleTab_.add( headRuleCharges_[i] );
        headRuleTab_.add( headRuleMSLevels_[i] ); 
        headRuleTab_.add( fragRuleMandatoryField );  
        headRuleTab_.add( deleteButton );      
      }
   }
   
    for(int i=0; i< headFragments_.size(); i++)
    {
    	for(int j=0; j< 3; j++)
    	{
    		refreshTextfieldToVerifieFirstTime(i, j, 1);
    	}
    }
   
   
   fragRuleOK_ = new JLabel( new ImageIcon( redCrossPicture_ ));
   gridBayLayout_.setConstraints(fragRuleOK_, new GridBagConstraints (0, headFragments_.size()+1, 1, 1, 0, 0, 
   GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0 ));
   headRuleTab_.add( fragRuleOK_ ); 
   fragRuleOK_.setToolTipText("Shows if the fragment is found in the spectra");
   
   fragRuleNameField_ = new FragmentNameTextField(null,new AddHeadFragmentDocumentListener(this));
   gridBayLayout_.setConstraints(fragRuleNameField_, new GridBagConstraints (1, headFragments_.size()+1, 1, 1, 0, 0, 
   GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0 ));
   fragRuleNameField_.setToolTipText("Enter the name of the fragment here!");   
   
   fragRuleFormulaField_ = new FragmentFormulaTextField(null,new AddHeadFragmentDocumentListener(this));
   gridBayLayout_.setConstraints(fragRuleFormulaField_, new GridBagConstraints (2, headFragments_.size()+1, 1, 1, 0, 0, 
   GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0 ));
   fragRuleFormulaField_.setToolTipText("Enter the formula of the fragment here!");
    
   fragRuleChargeField_ = new JTextField();
   gridBayLayout_.setConstraints(fragRuleChargeField_, new GridBagConstraints (3, headFragments_.size()+1, 1, 1, 0, 0, 
   GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0 ));
   fragRuleChargeField_.setColumns( 5 );
   fragRuleChargeField_.setToolTipText("Enter the charge of the fragment here!");
   fragRuleChargeField_.getDocument().addDocumentListener(new AddHeadFragmentDocumentListener(this));
    
   fragRuleMsLevelField_ = new JTextField();
   fragRuleMsLevelField_.setDocument(new IntegerRangeDocument(MIN_MSLEVEL,highestMSLevel_));
   gridBayLayout_.setConstraints(fragRuleMsLevelField_, new GridBagConstraints (4, headFragments_.size()+1, 1, 1, 0, 0, 
   GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0 ));
   fragRuleMsLevelField_.setColumns( 5 );   
   fragRuleMsLevelField_.setToolTipText("Enter the ms lever of the fragment here!");
   fragRuleMsLevelField_.getDocument().addDocumentListener(new AddHeadFragmentDocumentListener(this));
   
   fragRuleMandatoryCombo_ = new JComboBox<Object>(possibibilityListPostProcessing_); 
   fragRuleMandatoryCombo_.setSelectedIndex(1);   
   gridBayLayout_.setConstraints(fragRuleMandatoryCombo_, new GridBagConstraints (5, headFragments_.size()+1, 1, 1, 0, 0, 
   GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0 ));
   fragRuleMandatoryCombo_.setToolTipText("Choose the mandatory here!");
   fragRuleMandatoryCombo_.addActionListener(new ActionListener() 
   {
     public void actionPerformed(ActionEvent eventSource) 
     {           
       try {
        checkToAddHeadFragment();
      }
      catch (IOException | RulesException | SpectrummillParserException | NoRuleException | CgException |
          HydroxylationEncodingException | ChemicalFormulaException | LipidCombinameEncodingException e) {
        e.printStackTrace();
      }
     }
   });           

   
   JButton addFragmentFragRuleButton = new JButton( "Add Fragment" );
   gridBayLayout_.setConstraints(addFragmentFragRuleButton, new GridBagConstraints (1, headFragments_.size()+2, 1, 1, 0, 0, 
   GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0 ));
   addFragmentFragRuleButton.setToolTipText("Adds the fragment to the others");
   addFragmentFragRuleButton.addActionListener(new ActionListener() 
   {
     public void actionPerformed(ActionEvent event) 
     {       
      try 
      {
        FragmentRuleVO ruleVO = new FragmentRuleVO(fragRuleNameField_.getText(),fragRuleFormulaField_.getText(),Integer.parseInt(fragRuleChargeField_.getText()), 
        Integer.parseInt(fragRuleMsLevelField_.getText()), getFragMandatory(fragRuleMandatoryCombo_.getSelectedIndex()), headFragmentRules_, chainFragmentRules_, elementParser_);
        headFragmentRules_.put(fragRuleNameField_.getText(), ruleVO);        
        refreshMiddle(1);       
      }
      catch (RulesException e2) 
      {   
        JOptionPane.showMessageDialog(topSplitPane_, e2.getMessage());         
      }
      catch (java.lang.NumberFormatException e2) 
      {   
        JOptionPane.showMessageDialog(topSplitPane_, e2.getMessage());         
      }
     }
   }); 
 
   JLabel space1 = new JLabel( "     " );        
   gridBayLayout_.setConstraints(space1, new GridBagConstraints (1, headFragments_.size()+3, 1, 1, 0, 0, 
   GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0 ));
    
   JLabel equRuleEquation = new JLabel( "Equation" );        
   gridBayLayout_.setConstraints(equRuleEquation, new GridBagConstraints (1, headFragments_.size()+4, 1, 1, 0, 0, 
   GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0 ));
    
   JLabel equRuleMandatory = new JLabel( "Mandatory" );        
   gridBayLayout_.setConstraints(equRuleMandatory, new GridBagConstraints (2, headFragments_.size()+4, 1, 1, 0, 0, 
   GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0 ));
   
   int headIntensityRulesSize = headIntensityRules_.size();
   if(headIntensityRulesSize != 0)
   {
  	 headEquationVerifier_ = new EquationVerifier [headIntensityRules_.size()];
  	 headRuleEquationFields_ = new JTextField [headIntensityRules_.size()];
  	 headEquationMandatories_ = new boolean [headIntensityRules_.size()];

   for(int i=0; i < headIntensityRulesSize; i++)
   {    
  	 final int equationTablePosition = i; 
  	 
     EquationVerifier headEquationVerifier = new EquationVerifier(this, 1, equationTablePosition);     
     JTextField equRuleEquationField = headEquationVerifier.getTextField();
     equRuleEquationField.setColumns( 20 );
     gridBayLayout_.setConstraints(equRuleEquationField, new GridBagConstraints (1, headFragments_.size()+5+i, 1, 1, 0, 0, 
     GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0 ));     
     equRuleEquationField.setText( headIntensityRules_.elementAt(i).getRuleIdentifier() );
     equRuleEquationField.setEditable( true );      
     equRuleEquationField.getDocument().addDocumentListener(new ChangeEquationDocumentListener(this, equationTablePosition, 1));     
     headEquationVerifier_[i] = headEquationVerifier;
     headRuleEquationFields_[i] = equRuleEquationField;
      
  
     JComboBox<Object> equRuleMandatoryField = new JComboBox<Object>(possibibilityListPostProcessing_); 
     if(headIntensityRules_.elementAt(i).isMandatory())
     {
       equRuleMandatoryField.setSelectedIndex(0); 
     }           
     else
     {
       equRuleMandatoryField.setSelectedIndex(1);
     }        
     gridBayLayout_.setConstraints(equRuleMandatoryField, new GridBagConstraints (2, headFragments_.size()+5+i, 1, 1, 0, 0, 
     GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0 ));
     equRuleMandatoryField.addActionListener (new ActionListener () 
     {
       public void actionPerformed(ActionEvent e) 
       {
      	 headEquationMandatories_[equationTablePosition] = intToBoolean(equRuleMandatoryField.getSelectedIndex());
      	 refreshEquationFieldToVerifie(equationTablePosition, 1);       	       
       }
     });     
     headEquationMandatories_[equationTablePosition] = intToBoolean(equRuleMandatoryField.getSelectedIndex());  	  	     
     
     JButton deleteButton = new JButton();
     gridBayLayout_.setConstraints(deleteButton, new GridBagConstraints (3, headFragments_.size()+5+i, 1, 1, 0, 0, 
     GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0 ));
     deleteButton.setToolTipText("Deletes this equation");
     try
     {
     Image img = ImageIO.read(getClass().getResource("/images/error.gif"));
     deleteButton.setIcon(new ImageIcon(img));
     }
     catch (IOException e)
     {      
     }     
     final int currentElement = i;
     deleteButton.addActionListener(new ActionListener() 
     {
       public void actionPerformed(ActionEvent event) 
       {
         int reply = JOptionPane.showConfirmDialog(topSplitPane_, "Delete \"" + headIntensityRules_.elementAt(currentElement).getRuleIdentifier() 
             + "\" rule?", "Delete Rule" , JOptionPane.YES_NO_OPTION);
         if (reply == JOptionPane.YES_OPTION) 
         {
           headIntensityRules_.remove(currentElement);
           refreshMiddle(1);
           paintNewSpectra(false);  
         }                
       }
     });     
     headRuleTab_.add( equRuleEquationField );
     headRuleTab_.add( equRuleMandatoryField );
     headRuleTab_.add( deleteButton );
   }
   }
   
   //Refresh the Verifiers
   for(int i=0; i< headIntensityRules_.size(); i++)
   {
   		refreshEquationTextfieldToVerifieFirstTime(i,1);   	
   }
   
   equRuleOK_ =new JLabel( new ImageIcon( redCrossPicture_ ));
   gridBayLayout_.setConstraints(equRuleOK_, new GridBagConstraints (0, headIntensityRulesSize+headFragments_.size()+5, 1, 1, 0, 0, 
   GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0 ));
   equRuleOK_.setToolTipText("Shows if the equation is fulfilled");
    
   equRuleEquationField_ = new JTextField();
   gridBayLayout_.setConstraints(equRuleEquationField_, new GridBagConstraints (1, headIntensityRulesSize+headFragments_.size()+5, 1, 1, 0, 0, 
   GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0 ));
   equRuleEquationField_.setColumns( 20 );   
   equRuleEquationField_.setToolTipText("Enter the equation here!");
   equRuleEquationField_.getDocument().addDocumentListener(new AddHeadEquationDocumentListener(this));
      
   equRuleMandatoryCombo_ = new JComboBox<Object>(possibibilityListPostProcessing_); 
   equRuleMandatoryCombo_.setSelectedIndex(1);   
   gridBayLayout_.setConstraints(equRuleMandatoryCombo_, new GridBagConstraints (2, headIntensityRulesSize+headFragments_.size()+5, 1, 1, 0, 0,  
   GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0 ));   
   equRuleMandatoryCombo_.setToolTipText("Choose the mandatory here!");
   equRuleMandatoryCombo_.addActionListener (new ActionListener () 
   {
     public void actionPerformed(ActionEvent e) 
     {
       try {
        checkToAddHeadIntensityRule();
       } catch (RulesException e0){
         new WarningMessage(new JFrame(), "Error", e0.getMessage());
       } catch (IOException | SpectrummillParserException
           | NoRuleException | CgException | ChemicalFormulaException e1) {
         e1.printStackTrace();
       }             
     }
   });     

   
   JButton addFragmentEquRuleButton = new JButton( "Add Equation" );
   gridBayLayout_.setConstraints(addFragmentEquRuleButton, new GridBagConstraints (1, headIntensityRulesSize+headFragments_.size()+6, 1, 1, 0, 0, 
   GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0 )); 
   addFragmentEquRuleButton.setToolTipText("Adds the equation to the others");
   addFragmentEquRuleButton.addActionListener(new ActionListener() 
   {     
     public void actionPerformed(ActionEvent event) 
     {      
      try 
      {  
        IntensityRuleVO ruleVO = FragRuleParser.extractIntensityVOFromEquation(equRuleEquationField_.getText(), 0, FragRuleParser.HEAD_SECTION, FragmentRuleVO.getStringKeyHash(headFragmentRules_), FragmentRuleVO.getStringKeyHash(chainFragmentRules_),generalSettingsVO_.getAmountOfChains());
        ruleVO.setMandatory( getIntensityMandatory(equRuleMandatoryCombo_.getSelectedIndex()));
        headIntensityRules_.add(ruleVO);        
        refreshMiddle(1);  
        paintNewSpectra(false);
      }
      catch (RulesException e2) 
      {   
        JOptionPane.showMessageDialog(topSplitPane_, e2.getMessage());         
      }
     }    
   }); 
   
   
    
   JLabel space2 = new JLabel( "     " );        
   gridBayLayout_.setConstraints(space2, new GridBagConstraints (0, headIntensityRulesSize+headFragments_.size()+7, 1, 1, 0, 0, 
   GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0 ));
    
   JButton addShowDetailsButton = new JButton( "Show Details" );
   gridBayLayout_.setConstraints(addShowDetailsButton, new GridBagConstraints (1, headIntensityRulesSize+headFragments_.size()+8, 1, 1, 0, 0, 
   GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0 ));
   addShowDetailsButton.setToolTipText("Shows the details about the decision");
   addShowDetailsButton.addActionListener(new ActionListener() 
   {
     public void actionPerformed(ActionEvent event) 
     {             
       try 
       {     
         msnAnalyzer_ = updateMSnAnalyzerToCurrentSettings();
         Hashtable<String,CgProbe> detectedFragments = msnAnalyzer_.getHeadGroupFragments(); 
         Hashtable<String,IntensityRuleVO> fulfilledIntensityRules = msnAnalyzer_.getFulfilledHeadIntensityRules();
         printInDialogShowDetailsBox(detectedFragments, fulfilledIntensityRules, null, null);          
       } catch (RulesException rx) {       
         if (rx.getMessage().endsWith(FragRuleParser.NO_HEAD_AND_CHAINS_SECTION)){
           new WarningMessage(new JFrame(), "Error", "There are no rules entered!");
         } else rx.printStackTrace();;
       }       
       catch (IOException e) 
       {       
         e.printStackTrace();
       }       
       catch (SpectrummillParserException e) 
       {       
         e.printStackTrace();
       }
       catch (CgException e) 
       {       
         e.printStackTrace();
       }
       catch (NoRuleException e) 
       {
         e.printStackTrace();
       }
      catch (HydroxylationEncodingException | ChemicalFormulaException | LipidCombinameEncodingException e ) {
        e.printStackTrace();
      }      

     }     
   }); 
     
   headRuleTab_.add( fragRuleName );
   headRuleTab_.add( fragRuleFormula );
   headRuleTab_.add( fragRuleCharge );
   headRuleTab_.add( fragRuleMsLevel );
   headRuleTab_.add( fragRuleMandatory );
   headRuleTab_.add( fragRuleNameField_ );
   headRuleTab_.add( fragRuleFormulaField_ );
   headRuleTab_.add( fragRuleChargeField_ );
   headRuleTab_.add( fragRuleMsLevelField_ );
   headRuleTab_.add( fragRuleMandatoryCombo_ );
   headRuleTab_.add( addFragmentFragRuleButton );   
   headRuleTab_.add( space1 );
   headRuleTab_.add( equRuleOK_ );
   headRuleTab_.add( equRuleEquation );
   headRuleTab_.add( equRuleMandatory );
   headRuleTab_.add( equRuleEquationField_ );
   headRuleTab_.add( equRuleMandatoryCombo_ );
   headRuleTab_.add( addFragmentEquRuleButton );
   headRuleTab_.add( space2 );
   headRuleTab_.add( addShowDetailsButton );   
   headRuleTab_.setVisible(showFragmentTabs_);
   
   JScrollPane headRule = new JScrollPane(headRuleTab_); 
   middleTabSection_.addTab("Head Rules", headRule);
   
  }
  
  /**
   * Delets the possible dummy fragment
   */
  private void deletePossibleDummy()
	{
  	if(headFragmentRules_.containsKey("EMPTY"))
  	{
			headFragmentRules_.remove("EMPTY"); 
  	}
  	headFragments_ = getSortedVector(headFragmentRules_, TYPE_HEAD);	
	}

	/**
   * Creates the chain rule tab and adds it to the middle tab section 
   */
  public void createChainTab()
  {   
    chainTab_ = new JPanel();   
    chainTab_.setLayout(gridBayLayout_);  
 
    JLabel chainRuleName = new JLabel( "Name" );        
    gridBayLayout_.setConstraints(chainRuleName, new GridBagConstraints (1, 0, 1, 1, 0, 0, 
    GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0 ));
    
    JLabel chainRuleFormula = new JLabel( "Formula" );        
    gridBayLayout_.setConstraints(chainRuleFormula, new GridBagConstraints (2, 0, 1, 1, 0, 0, 
    GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0 ));
    
    JLabel chainRuleCharge = new JLabel( "Charge" );        
    gridBayLayout_.setConstraints(chainRuleCharge, new GridBagConstraints (3, 0, 1, 1, 0, 0, 
    GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0 ));
    
    JLabel chainRuleMsLevel = new JLabel( "MS Level" );        
    gridBayLayout_.setConstraints(chainRuleMsLevel, new GridBagConstraints (4, 0, 1, 1, 0, 0, 
    GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0 ));
    
    JLabel chainRuleMandatory = new JLabel( "Mandatory" );        
    gridBayLayout_.setConstraints(chainRuleMandatory, new GridBagConstraints (5, 0, 1, 1, 0, 0, 
    GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0 ));
      
    sortedChainFragments_ = getSortedVector(chainFragmentRules_, TYPE_CHAIN);
    if(sortedChainFragments_.size() !=0)
    {
      chainFragmentFormulaVerifier_ = new FragmentVerifier [chainFragmentRules_.size()];
      chainRuleFormulas_ = new JTextField [chainFragmentRules_.size()];
      
      chainFragmentChargeVerifier_ = new FragmentVerifier [chainFragmentRules_.size()];
      chainRuleCharges_ = new JTextField [chainFragmentRules_.size()];
      
      chainFragmentMsLevelVerifier_ = new FragmentVerifier [chainFragmentRules_.size()];
      chainRuleMSLevels_ = new JTextField [chainFragmentRules_.size()];
    
      chainFragmentMandatories_ = new boolean [chainFragmentRules_.size()];
    
    for(int i=0; i<sortedChainFragments_.size(); i++ )
    {     	
    	final int positionInTable = i; 
    	
    	chainFragRuleName_ = new JTextField();
      gridBayLayout_.setConstraints(chainFragRuleName_, new GridBagConstraints (1, i+1, 1, 1, 0, 0, 
      GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0 ));
      chainFragRuleName_.setColumns( 20 );
      chainFragRuleName_.setText( sortedChainFragments_.elementAt(i).getName() );
      chainFragRuleName_.setEditable(false);
      
      FragmentVerifier chainFragmentFormulaVerifier = new FragmentVerifier(this, FragmentVerifier.TYPE_FORMULA, 2);
      JTextField chainRuleFormulaField = chainFragmentFormulaVerifier.getTextField();
      gridBayLayout_.setConstraints(chainRuleFormulaField, new GridBagConstraints (2, i+1, 1, 1, 0, 0, 
      GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0 ));      
      chainRuleFormulaField.setText( sortedChainFragments_.elementAt(i).getFormula() );
      chainRuleFormulaField.getDocument().addDocumentListener(new ChangeFragmentDocumentListener(this, positionInTable, 0, 2)); 
      chainRuleFormulaField.addFocusListener(new ChangeFragmentFocusListener(this, positionInTable, 0, 2));     
      chainFragmentFormulaVerifier_[i] = chainFragmentFormulaVerifier;
      chainRuleFormulas_[i] = chainRuleFormulaField;
       
      FragmentVerifier chainFragmentChargeVerifier = new FragmentVerifier(this, 2, 2);
      JTextField chainRuleChargeField = chainFragmentChargeVerifier.getTextField();
      chainRuleChargeField.setColumns( 5 );
      gridBayLayout_.setConstraints(chainRuleChargeField, new GridBagConstraints (3, i+1, 1, 1, 0, 0, 
      GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0 ));      
      chainRuleChargeField.setText( Integer.toString(sortedChainFragments_.elementAt(i).getCharge()) );
      chainRuleChargeField.setEditable(true);      
      chainRuleChargeField.getDocument().addDocumentListener(new ChangeFragmentDocumentListener(this, positionInTable, 1, 2));     
      chainRuleChargeField.addFocusListener(new ChangeFragmentFocusListener(this, positionInTable, 1, 2));     
      chainFragmentChargeVerifier_[i] = chainFragmentChargeVerifier;
      chainRuleCharges_[i] = chainRuleChargeField;      
       
      FragmentVerifier chainFragmentMsLevelVerifier = new FragmentVerifier(this, 3, 2);
      JTextField chainRuleMsLevelField = chainFragmentMsLevelVerifier.getTextField();
      chainRuleMsLevelField.setDocument(new IntegerRangeDocument(MIN_MSLEVEL,highestMSLevel_));
      chainRuleMsLevelField.setColumns( 5 );
      gridBayLayout_.setConstraints(chainRuleMsLevelField, new GridBagConstraints (4, i+1, 1, 1, 0, 0, 
      GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0 ));      
      chainRuleMsLevelField.setText( Integer.toString(sortedChainFragments_.elementAt(i).getMsLevel()) );
      chainRuleMsLevelField.setEditable(true);      
      chainRuleMsLevelField.getDocument().addDocumentListener(new ChangeFragmentDocumentListener(this, positionInTable, 2, 2));
      chainRuleMsLevelField.addFocusListener(new ChangeFragmentFocusListener(this, positionInTable, 2, 2)); 
      chainFragmentMsLevelVerifier_[i] = chainFragmentMsLevelVerifier;
      chainRuleMSLevels_[i] = chainRuleMsLevelField;
           
      final JComboBox<Object> chainRuleMandatoryField = new JComboBox<Object>(possibibilityListPostProcessing_); 
      if(sortedChainFragments_.elementAt(i).isMandatory()==FragmentRuleVO.MANDATORY_TRUE)
      {
        chainRuleMandatoryField.setSelectedIndex(0); 
      }           
      else
      {
        chainRuleMandatoryField.setSelectedIndex(1);
      }        
      gridBayLayout_.setConstraints(chainRuleMandatoryField, new GridBagConstraints (5, i+1, 1, 1, 0, 0, 
      GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0 ));
      chainRuleMandatoryField.addActionListener (new ActionListener () 
      {
      	public void actionPerformed(ActionEvent e) 
      	{        	
      		chainFragmentMandatories_[positionInTable] = intToBoolean(chainRuleMandatoryField.getSelectedIndex());
        	refreshTextfieldToVerifie(positionInTable, 0, 2);
        	checkFragmentWithErrors(positionInTable, 2, 2);            
      	}
      });     
       
      JButton deleteButton = new JButton();
      gridBayLayout_.setConstraints(deleteButton, new GridBagConstraints (6, i+1, 1, 1, 0, 0, 
      GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0 ));  
      deleteButton.setToolTipText("Deletes this fragment");       
      
      try
      {
        Image img = ImageIO.read(getClass().getResource("/images/error.gif"));
        deleteButton.setIcon(new ImageIcon(img));
      }
      catch (IOException e1)
      {      
      }     
      final String key = sortedChainFragments_.elementAt(i).getName();
      deleteButton.addActionListener(new ActionListener() 
      {
        public void actionPerformed(ActionEvent event) 
        {
          int reply = JOptionPane.showConfirmDialog(topSplitPane_, "Delete \"" + chainFragmentRules_.get(key).getName() + "\" rule?", "Delete Rule" , JOptionPane.YES_NO_OPTION);
          if (reply == JOptionPane.YES_OPTION) 
          {
            chainFragmentRules_.remove(key);
            deleteDependingRules(key,TYPE_CHAIN);
            refreshMiddle(2);
            paintNewSpectra(false);  
          }                
        }
      });
       
      chainTab_.add(chainFragRuleName_);
      chainTab_.add(chainRuleFormulas_[i]);
      chainTab_.add(chainRuleCharges_[i]);
      chainTab_.add(chainRuleMSLevels_[i]);     
      chainTab_.add(chainRuleMandatoryField);
      chainTab_.add(deleteButton);
    }
    }
    
    for(int i=0; i< sortedChainFragments_.size(); i++)
    {
    	for(int j=0; j< 3; j++)
    	{
    		refreshTextfieldToVerifieFirstTime(i, j, 2);
    	}
    }
    
    
    chainRuleOK_ = new JLabel( new ImageIcon( redCrossPicture_ ));    
    gridBayLayout_.setConstraints(chainRuleOK_, new GridBagConstraints (0, chainFragmentRules_.size()+1, 1, 1, 0, 0, 
    GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0 ));
    chainTab_.add(chainRuleOK_);
    chainRuleOK_.setToolTipText("Shows if the fragment is found in the spectra");
     
    chainRuleNameField_ = new FragmentNameTextField(null,new AddChainFragmentDocumentListener(this));
    gridBayLayout_.setConstraints(chainRuleNameField_, new GridBagConstraints (1, chainFragmentRules_.size()+1, 1, 1, 0, 0, 
    GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0 ));
    chainRuleNameField_.setToolTipText("Enter name of the fragment here!");
    
     
    chainRuleFormulaField_ = new FragmentFormulaTextField(null,new AddChainFragmentDocumentListener(this));
    gridBayLayout_.setConstraints(chainRuleFormulaField_, new GridBagConstraints (2, chainFragmentRules_.size()+1, 1, 1, 0, 0, 
    GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0 ));
    chainRuleFormulaField_.setToolTipText("Enter formula of the fragment here!");
    
     
    chainRuleChargeField_ = new JTextField();
    gridBayLayout_.setConstraints(chainRuleChargeField_, new GridBagConstraints (3, chainFragmentRules_.size()+1, 1, 1, 0, 0, 
    GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0 ));
    chainRuleChargeField_.setColumns( 5 );
    chainRuleChargeField_.setToolTipText("Enter the charge of the fragment here!");
    chainRuleChargeField_.getDocument().addDocumentListener(new AddChainFragmentDocumentListener(this));
     
    chainRuleMsLevelField_ = new JTextField();
    chainRuleMsLevelField_.setDocument(new IntegerRangeDocument(MIN_MSLEVEL,highestMSLevel_));
    gridBayLayout_.setConstraints(chainRuleMsLevelField_, new GridBagConstraints (4, chainFragmentRules_.size()+1, 1, 1, 0, 0, 
    GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0 ));
    chainRuleMsLevelField_.setColumns( 5 );
    chainRuleMsLevelField_.setVisible(true);
    chainRuleMsLevelField_.setToolTipText("Enter the ms lever of the fragment here!");
    chainRuleMsLevelField_.getDocument().addDocumentListener(new AddChainFragmentDocumentListener(this));

    chainRuleMandatoryCombo_ = new JComboBox<Object>(possibibilityListPostProcessing_);
    chainRuleMandatoryCombo_.setSelectedIndex(1); 
    gridBayLayout_.setConstraints(chainRuleMandatoryCombo_, new GridBagConstraints (5, chainFragmentRules_.size()+1, 1, 1, 0, 0, 
    GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0 ));
    chainRuleMandatoryCombo_.setToolTipText("Choose the mandatory here!"); 
    chainRuleMandatoryCombo_.addActionListener (new ActionListener () 
    {
    public void actionPerformed(ActionEvent e) 
    {        
      try {
        checkToAddChainFragment();
      }
      catch (RulesException | IOException | SpectrummillParserException
          | CgException | NoRuleException | HydroxylationEncodingException | ChemicalFormulaException | LipidCombinameEncodingException e1) {
        e1.printStackTrace();
      }      
    }
    });
    
    JButton addFragmentChainRuleButton = new JButton( "Add Fragment" );
    gridBayLayout_.setConstraints(addFragmentChainRuleButton, new GridBagConstraints (1, chainFragmentRules_.size()+2, 1, 1, 0, 0, 
    GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0 ));
    addFragmentChainRuleButton.setToolTipText("Adds the fragment to the others");
    addFragmentChainRuleButton.addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent event) 
      {
       
       try 
       {
         FragmentRuleVO ruleVO = new FragmentRuleVO(chainRuleNameField_.getText(),chainRuleFormulaField_.getText(),Integer.parseInt(chainRuleChargeField_.getText()), 
         Integer.parseInt(chainRuleMsLevelField_.getText()), getFragMandatory(chainRuleMandatoryCombo_.getSelectedIndex()), headFragmentRules_, chainFragmentRules_, elementParser_);
         chainFragmentRules_.put(chainRuleNameField_.getText(), ruleVO);
         refreshMiddle(2);
        
       }
       catch (RulesException e2) 
       {   
         JOptionPane.showMessageDialog(topSplitPane_, e2.getMessage());         
       }
       catch (java.lang.NumberFormatException e2) 
       {   
         JOptionPane.showMessageDialog(topSplitPane_, e2.getMessage());         
       }
      }
    });
     
    JLabel space4 = new JLabel( "     " );        
    gridBayLayout_.setConstraints(space4, new GridBagConstraints (1, chainFragmentRules_.size()+3, 1, 1, 0, 0, 
    GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0 ));   
           
    JLabel equChainRuleEquation = new JLabel( "Equation" );        
    gridBayLayout_.setConstraints(equChainRuleEquation, new GridBagConstraints (1, chainFragmentRules_.size()+4, 1, 1, 0, 0, 
    GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0 ));
      
    JLabel equChainRuleMandatory = new JLabel( "Mandatory" );        
    gridBayLayout_.setConstraints(equChainRuleMandatory, new GridBagConstraints (2, chainFragmentRules_.size()+4, 1, 1, 0, 0, 
    GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0 ));
    
    int chainIntensityRulesSize = chainIntensityRules_.size();
    if(chainIntensityRulesSize != 0)
    {
    	chainEquationVerifier_ = new EquationVerifier [chainIntensityRules_.size()];
   	 	chainRuleEquationFields_ = new JTextField [chainIntensityRules_.size()];
   	 	chainEquationMandatories_ = new boolean [chainIntensityRules_.size()];
   	 
    for(int i=0; i < chainIntensityRulesSize; i++)
    {
      final int equationTablePosition = i;
      
      
      EquationVerifier chainEquationVerifier = new EquationVerifier(this, 2, equationTablePosition);      
      JTextField equChainRuleEquationField = chainEquationVerifier.getTextField();
      equChainRuleEquationField.setColumns( 20 );
      gridBayLayout_.setConstraints(equChainRuleEquationField, new GridBagConstraints (1, chainFragmentRules_.size()+5+i, 1, 1, 0, 0, 
      GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0 ));      
      equChainRuleEquationField.setText( chainIntensityRules_.elementAt(i).getRuleIdentifier() );
      equChainRuleEquationField.setEditable(true);       
      equChainRuleEquationField.getDocument().addDocumentListener(new ChangeEquationDocumentListener(this, equationTablePosition, 2));
      chainEquationVerifier_[i] = chainEquationVerifier;
      chainRuleEquationFields_[i] = equChainRuleEquationField;
          
      
      JComboBox<Object> equChainRuleMandatoryField = new JComboBox<Object>(possibibilityListPostProcessing_); 
      if(chainIntensityRules_.elementAt(i).isMandatory())
      {
        equChainRuleMandatoryField.setSelectedIndex(0); 
      }           
      else
      {
        equChainRuleMandatoryField.setSelectedIndex(1);
      }        
      gridBayLayout_.setConstraints(equChainRuleMandatoryField, new GridBagConstraints (2, chainFragmentRules_.size()+5+i, 1, 1, 0, 0, 
          GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0 ));
      equChainRuleMandatoryField.addActionListener (new ActionListener () 
      {
        public void actionPerformed(ActionEvent e) 
        {
        	chainEquationMandatories_[equationTablePosition] = intToBoolean(equChainRuleMandatoryField.getSelectedIndex());       	 	
       	 	refreshEquationFieldToVerifie(equationTablePosition, 2);   
        }
      });
      chainEquationMandatories_[equationTablePosition] = intToBoolean(equChainRuleMandatoryField.getSelectedIndex());   	 	
      
      JButton deleteButton = new JButton();
      gridBayLayout_.setConstraints(deleteButton, new GridBagConstraints (3, chainFragmentRules_.size()+5+i, 1, 1, 0, 0, 
      GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0 ));
      deleteButton.setToolTipText("Deletes this equation");
      try
      {
      Image img = ImageIO.read(getClass().getResource("/images/error.gif"));
      deleteButton.setIcon(new ImageIcon(img));
      }
      catch (IOException e1)
      {      
      }
      final int currentElement = i;
      deleteButton.addActionListener(new ActionListener() 
      {
        public void actionPerformed(ActionEvent event) 
        {
          int reply = JOptionPane.showConfirmDialog(topSplitPane_, "Delete \"" + chainIntensityRules_.elementAt(currentElement).getRuleIdentifier() + "\" rule?", "Delete Rule" , JOptionPane.YES_NO_OPTION);
          if (reply == JOptionPane.YES_OPTION) 
          {
            chainIntensityRules_.remove(currentElement);
            refreshMiddle(2);
            paintNewSpectra(false);      
          }                
        }
      });
      
      
      chainTab_.add( equChainRuleEquationField );
      chainTab_.add( equChainRuleMandatoryField );
      chainTab_.add( deleteButton );
    }
    }
    
    equChainRuleOK_ =new JLabel( new ImageIcon( redCrossPicture_ ));
    gridBayLayout_.setConstraints(equChainRuleOK_, new GridBagConstraints (0, chainIntensityRulesSize+chainFragmentRules_.size()+5, 1, 1, 0, 0, 
    GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0 ));
    equChainRuleOK_.setToolTipText("Shows if the equation is fulfilled");
      
    equChainRuleEquationField_ = new JTextField();
    gridBayLayout_.setConstraints(equChainRuleEquationField_, new GridBagConstraints (1, chainIntensityRulesSize+chainFragmentRules_.size()+5, 1, 1, 0, 0, 
    GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0 ));
    equChainRuleEquationField_.setColumns( 20 );
    equChainRuleEquationField_.setToolTipText("Enter the equation here!");
    equChainRuleEquationField_.getDocument().addDocumentListener(new AddChainEquationDocumentListener(this)); 
       
    
    equChainRuleMandatoryCombo_ = new JComboBox<Object>(possibibilityListPostProcessing_); 
    equChainRuleMandatoryCombo_.setSelectedIndex(1);   
    gridBayLayout_.setConstraints(equChainRuleMandatoryCombo_, new GridBagConstraints (2, chainIntensityRulesSize+chainFragmentRules_.size()+5, 1, 1, 0, 0, 
    GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0 ));    
    equChainRuleMandatoryCombo_.setToolTipText("Choose the mandatory here!");
    equChainRuleMandatoryCombo_.addActionListener (new ActionListener () 
    {
      public void actionPerformed(ActionEvent e) 
      {
        try {
          checkToAddChainIntensityRule();
        } catch (RulesException e0){
          new WarningMessage(new JFrame(), "Error", e0.getMessage());
        } catch (IOException | SpectrummillParserException
            | NoRuleException | CgException | ChemicalFormulaException e1 ) {
          e1.printStackTrace();
        }
      }
    });

    
    JButton addFragmentChainRuleButton2 = new JButton( "Add Equation" );
    gridBayLayout_.setConstraints(addFragmentChainRuleButton2, new GridBagConstraints (1, chainIntensityRulesSize+chainFragmentRules_.size()+6, 1, 1, 0, 0, 
    GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0 ));
    addFragmentChainRuleButton2.setToolTipText("Adds the equation to the others");
    addFragmentChainRuleButton2.addActionListener(new ActionListener() 
    {     
      public void actionPerformed(ActionEvent event) 
      { 
           
       try 
       {
         IntensityRuleVO ruleVO = FragRuleParser.extractIntensityVOFromEquation(equChainRuleEquationField_.getText(), 0, FragRuleParser.CHAINS_SECTION, FragmentRuleVO.getStringKeyHash(headFragmentRules_), FragmentRuleVO.getStringKeyHash(chainFragmentRules_),generalSettingsVO_.getAmountOfChains());
         ruleVO.setMandatory(getIntensityMandatory(equChainRuleMandatoryCombo_.getSelectedIndex()));
         chainIntensityRules_.add(ruleVO);        
         refreshMiddle(2);    
         paintNewSpectra(false);   
       }
       catch (RulesException e2) 
       {   
         JOptionPane.showMessageDialog(topSplitPane_, e2.getMessage());         
       }
      }
     
    }); 
     
    JLabel space5 = new JLabel( "     " );        
    gridBayLayout_.setConstraints(space5, new GridBagConstraints (0, chainIntensityRulesSize+chainFragmentRules_.size()+7, 1, 1, 0, 0, 
    GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0 ));
    
    JButton addChainShowDetailsButton = new JButton( "Show Details" );
    gridBayLayout_.setConstraints(addChainShowDetailsButton, new GridBagConstraints (1, chainIntensityRulesSize+chainFragmentRules_.size()+8, 1, 1, 0, 0, 
    GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0 ));
    addChainShowDetailsButton.setToolTipText("Shows the details about the decision");
    addChainShowDetailsButton.addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent event) 
      {           
        try 
        {
          try 
          {
        	    msnAnalyzer_ = updateMSnAnalyzerToCurrentSettings();
          }
          catch (NoRuleException | HydroxylationEncodingException | ChemicalFormulaException | LipidCombinameEncodingException e) 
          {
            e.printStackTrace();
          }         
          Hashtable<String,Hashtable<String,CgProbe>> chainFragments = msnAnalyzer_.getChainFragments();           
          Hashtable<String,Hashtable<String,IntensityChainVO>> fulfilledChainIntensityRules = msnAnalyzer_.getFulfilledChainIntensityRules();
          printInDialogShowDetailsBox(null, null, chainFragments, fulfilledChainIntensityRules);          
        } catch (RulesException rx) {       
          if (rx.getMessage().endsWith(FragRuleParser.NO_HEAD_AND_CHAINS_SECTION)){
            new WarningMessage(new JFrame(), "Error", "There are no rules entered!");
          } else rx.printStackTrace();;
        }       
        catch (IOException e) 
        {       
          e.printStackTrace();
        }       
        catch (SpectrummillParserException e) 
        {       
          e.printStackTrace();
        }
        catch (CgException | ChemicalFormulaException | HydroxylationEncodingException | LipidCombinameEncodingException e) 
        {       
          e.printStackTrace();
        }       
      }     
    }); 
    
    chainTab_.add(chainRuleName);
    chainTab_.add(chainRuleFormula);
    chainTab_.add(chainRuleCharge);
    chainTab_.add(chainRuleMsLevel);
    chainTab_.add(chainRuleMandatory);    
    chainTab_.add(chainRuleNameField_);
    chainTab_.add(chainRuleFormulaField_);
    chainTab_.add(chainRuleChargeField_);
    chainTab_.add(chainRuleMsLevelField_);     
    chainTab_.add(chainRuleMandatoryCombo_);
    chainTab_.add(addFragmentChainRuleButton);
    chainTab_.add(space4);     
    chainTab_.add(equChainRuleEquation );
    chainTab_.add(equChainRuleMandatory);
    chainTab_.add(equChainRuleOK_);    
    chainTab_.add(equChainRuleEquationField_ );
    chainTab_.add(equChainRuleMandatoryCombo_ );
    chainTab_.add(addFragmentChainRuleButton2);
    chainTab_.add(space5);     
    chainTab_.add(addChainShowDetailsButton);
    chainTab_.setVisible(showFragmentTabs_);
    JScrollPane chain = new JScrollPane(chainTab_);    
    middleTabSection_.addTab("Chain", chain);  
  }
    
  /**
   * Creates the position tab and adds it to the middle tab section 
   */
  public void createPositionTab()
  {    
    positionTab_ = new JPanel();
    positionTab_.setLayout(gridBayLayout_);
      
    JLabel positionRuleEquation = new JLabel( "Equation" );        
    gridBayLayout_.setConstraints(positionRuleEquation, new GridBagConstraints (1, 0, 1, 1, 0, 0, 
        GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0 ));
    
    JLabel positionRuleMandatory = new JLabel( "Mandatory" );        
    gridBayLayout_.setConstraints(positionRuleMandatory, new GridBagConstraints (2, 0, 1, 1, 0, 0, 
        GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0 ));
    
 
    int positionIntensityRulesSize = positionIntensityRules_.size();    
    if(positionIntensityRulesSize != 0)
    {
    	
    	positionEquationVerifier_ = new EquationVerifier [positionIntensityRules_.size()];
   	 	positionRuleEquationFields_ = new JTextField [positionIntensityRules_.size()];
   	 	positionEquationMandatories_ = new boolean [positionIntensityRules_.size()];
   	 	
    for(int i=0; i < positionIntensityRulesSize; i++)
    {       
      final int equationTablePosition = i;
      
      EquationVerifier positionEquationVerifier = new EquationVerifier(this, 3, equationTablePosition);
      JTextField positionEquationField = positionEquationVerifier.getTextField();
      positionEquationField.setColumns( 20 );
      gridBayLayout_.setConstraints(positionEquationField, new GridBagConstraints (1, 1+i, 1, 1, 0, 0, 
      GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0 ));      
      positionEquationField.setText( positionIntensityRules_.elementAt(i).getRuleIdentifier() );
      positionEquationField.setEditable(true);
      positionEquationField.getDocument().addDocumentListener(new ChangeEquationDocumentListener(this, equationTablePosition, 3));
      positionEquationVerifier_[i] = positionEquationVerifier;
      positionRuleEquationFields_[i] = positionEquationField;
    
      JComboBox<Object> positionMandatoryField = new JComboBox<Object>(possibibilityListPostProcessing_); 
      if(positionIntensityRules_.elementAt(i).isMandatory())
      {
        positionMandatoryField.setSelectedIndex(0); 
      }           
      else
      {
        positionMandatoryField.setSelectedIndex(1);
      }        
      gridBayLayout_.setConstraints(positionMandatoryField, new GridBagConstraints (2, 1+i, 1, 1, 0, 0, 
      GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0 ));
      positionMandatoryField.addActionListener (new ActionListener () 
      {
        public void actionPerformed(ActionEvent e) 
        {
        	positionEquationMandatories_[equationTablePosition] = intToBoolean(positionMandatoryField.getSelectedIndex());
       	 	refreshEquationFieldToVerifie(equationTablePosition, 3);   
        }
      });      
      positionEquationMandatories_[equationTablePosition] = intToBoolean(positionMandatoryField.getSelectedIndex());
      
      JButton deleteButton = new JButton();
      gridBayLayout_.setConstraints(deleteButton, new GridBagConstraints (3, 1+i, 1, 1, 0, 0, 
      GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0 ));
      deleteButton.setToolTipText("Deletes this equation");
      try
      {
        Image img = ImageIO.read(getClass().getResource("/images/error.gif"));
        deleteButton.setIcon(new ImageIcon(img));
      }
      catch (IOException e1)
      {      
      }
      final int currentElement = i;
      deleteButton.addActionListener(new ActionListener() 
      {
        public void actionPerformed(ActionEvent event) 
        {
          int reply = JOptionPane.showConfirmDialog(topSplitPane_, "Delete \"" + positionIntensityRules_.elementAt(currentElement).getRuleIdentifier() + "\" rule?", "Delete Rule" , JOptionPane.YES_NO_OPTION);
          if (reply == JOptionPane.YES_OPTION) 
          {
            positionIntensityRules_.remove(currentElement);
            refreshMiddle(3);           
            paintNewSpectra(false);           
          }                
        }
      });
      
      positionTab_.add( positionEquationField );
      positionTab_.add( positionMandatoryField );
      positionTab_.add( deleteButton );
    }
    }
    
    
    positionRuleOK_ = new JLabel( new ImageIcon( redCrossPicture_ ));
    gridBayLayout_.setConstraints(positionRuleOK_, new GridBagConstraints (0, 1+positionIntensityRulesSize, 1, 1, 0, 0, 
    GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0 ));
    positionRuleOK_.setToolTipText("Shows if the equation is fulfilled");
    
    positionEquationField_ = new JTextField();
    gridBayLayout_.setConstraints(positionEquationField_, new GridBagConstraints (1, 1+positionIntensityRulesSize, 1, 1, 0, 0, 
    GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0 ));
    positionEquationField_.setColumns( 20 );
    positionEquationField_.setToolTipText("Enter the equation here!");
    positionEquationField_.getDocument().addDocumentListener(new AddPositionEquationDocumentListener(this));
    
    positionRuleMandatoryCombo_ = new JComboBox<Object>(possibibilityListPostProcessing_); 
    positionRuleMandatoryCombo_.setSelectedIndex(1);   
    gridBayLayout_.setConstraints(positionRuleMandatoryCombo_, new GridBagConstraints (2, 1+positionIntensityRulesSize, 1, 1, 0, 0, 
    GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0 ));
    positionRuleMandatoryCombo_.setToolTipText("Choose the mandatory here!");
    positionRuleMandatoryCombo_.addActionListener (new ActionListener () 
    {
      public void actionPerformed(ActionEvent e) 
      {
        try {
          checkToAddPositionRule();
        }
        catch (RulesException | IOException | SpectrummillParserException
            | NoRuleException | CgException e1) {
          // TODO Auto-generated catch block
          e1.printStackTrace();
        }
      }
    });      

    
    JButton addFragmentPositionButton = new JButton( "Add Equation" );
    gridBayLayout_.setConstraints(addFragmentPositionButton, new GridBagConstraints (1, 6+positionIntensityRulesSize, 1, 1, 0, 0, 
    GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0 ));
    addFragmentPositionButton.setToolTipText("Adds the equation to the others");
    addFragmentPositionButton.addActionListener(new ActionListener() 
    {     
      public void actionPerformed(ActionEvent event) 
      {         
       try 
       {
         IntensityRuleVO ruleVO = FragRuleParser.extractIntensityVOFromEquation(positionEquationField_.getText(), 0, FragRuleParser.POSITION_SECTION, FragmentRuleVO.getStringKeyHash(headFragmentRules_), FragmentRuleVO.getStringKeyHash(chainFragmentRules_),generalSettingsVO_.getAmountOfChains());         
         ruleVO.setMandatory(getIntensityMandatory(positionRuleMandatoryCombo_.getSelectedIndex()));
         positionIntensityRules_.add(ruleVO);        
         refreshMiddle(3);   
         paintNewSpectra(false);
       }
       catch (RulesException e2) 
       {   
         JOptionPane.showMessageDialog(topSplitPane_, e2.getMessage());         
       } 
      }
     
    }); 
    
    JLabel space3 = new JLabel( "     " );        
    gridBayLayout_.setConstraints(space3, new GridBagConstraints (0, 7+positionIntensityRulesSize, 1, 1, 0, 0, 
    GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0 ));
    
    JButton addPositionShowDetailsButton = new JButton( "Show Details" );
    gridBayLayout_.setConstraints(addPositionShowDetailsButton, new GridBagConstraints (1, 8+positionIntensityRulesSize, 1, 1, 0, 0, 
    GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0 ));
    addPositionShowDetailsButton.setToolTipText("Shows the details about the decision");
    addPositionShowDetailsButton.addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent event) 
      {            
        try 
        {         
          msnAnalyzer_ = updateMSnAnalyzerToCurrentSettings();
          printPositionRecomandation(msnAnalyzer_.getResult(),msnAnalyzer_.getPositionDefinition());
        } catch (RulesException rx) {       
          if (rx.getMessage().endsWith(FragRuleParser.NO_HEAD_AND_CHAINS_SECTION)){
            new WarningMessage(new JFrame(), "Error", "There are no rules entered!");
          } else rx.printStackTrace();;
        }       
        catch (IOException e) 
        {       
          e.printStackTrace();
        }       
        catch (SpectrummillParserException e) 
        {       
          e.printStackTrace();
        }
        catch (CgException e) 
        {       
          e.printStackTrace();
        }
        catch (NoRuleException e) 
        {
          e.printStackTrace();
        }
        catch (HydroxylationEncodingException | ChemicalFormulaException | LipidCombinameEncodingException e) {
          e.printStackTrace();
        }  
      }     
    }); 
    
    
    
    positionTab_.add( positionRuleEquation );
    positionTab_.add( positionRuleMandatory );
    positionTab_.add( positionRuleOK_ );
    positionTab_.add( positionEquationField_ );
    positionTab_.add( positionRuleMandatoryCombo_ );
    positionTab_.add( addFragmentPositionButton );
    positionTab_.add( space3 );
    positionTab_.add( addPositionShowDetailsButton );    
    JScrollPane position = new JScrollPane(positionTab_);  
    positionTab_.setVisible(showFragmentTabs_);
    middleTabSection_.addTab("Position", position);    
  }

  
  
  /**
   * Creates the bottom section where the changes can be submitted and checked
   */
  public void finalButtonsSection() 
  {
    bottomSectionButtons_ = new JPanel(); 
    bottomSectionButtons_.setLayout(gridBayLayout_);    
   
    JButton showTotalDecisionButton = new JButton( "Show Total Decision" );
    gridBayLayout_.setConstraints(showTotalDecisionButton, new GridBagConstraints (0, 0, 1, 1, 0, 0, 
    GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0 )); 
    showTotalDecisionButton.setToolTipText("Shows the total deciosion");
    showTotalDecisionButton.addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent e)
      {  
        printTotalDecisionInBox();       
      }
    });    
    
    JButton saveRulesButton = new JButton( "Save Rules" );
    gridBayLayout_.setConstraints(saveRulesButton, new GridBagConstraints (1, 0, 1, 1, 0, 0, 
    GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0 ));    
    saveRulesButton.setVisible(true);
    saveRulesButton.setToolTipText("Saves the rules in the .txt file");
    saveRulesButton.addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent event) 
      {        
        int reply = JOptionPane.showConfirmDialog(topSplitPane_, "Do you want to save the changed Rules? Existing Rules are overwritten and influence the running calculations." , "Save Rules" , JOptionPane.YES_NO_OPTION);
        if (reply == JOptionPane.YES_OPTION) 
        {
          String fileBeforeSave;
          try 
          {
            fileBeforeSave = readFile("fragRules/"+lipidClassName_ + "_" + lipidAdduct_ + ".frag.txt");
            if (fileBeforeSave.length() > 0) 
            {
            	fileBeforeSave = fileBeforeSave.substring(0, fileBeforeSave.length()-1);
            }
          }
          catch (IOException e1) 
          {
            fileBeforeSave = "EMPTY";
          }
        
            try
            {      
            	if(generalSettings_.checkPercentageFieldsForNumbers())
            	{
              setRuleClassIdentifier(); 
              saveRules("default");             
	              getRules();
	              JOptionPane.showMessageDialog(topSplitPane_, "Rules successfully saved.");
	              refreshMiddle(0); 
            	}
            }
            catch (IOException e) 
            {
              JOptionPane.showMessageDialog(topSplitPane_, e.getMessage());
              try 
              {
                writeStringInCurrentFile(fileBeforeSave);
              }
              catch (IOException e1) 
              {
              }
            } 
            catch (SpectrummillParserException e) 
            {
              JOptionPane.showMessageDialog(topSplitPane_, e.getMessage());   
              try 
              {
                writeStringInCurrentFile(fileBeforeSave);
              }
              catch (IOException e1) 
              {
              }
            }
            catch (RulesException e) 
            {
              JOptionPane.showMessageDialog(topSplitPane_, e.getMessage());  
              try 
              {
                writeStringInCurrentFile(fileBeforeSave);
              }
              catch (IOException e1) 
              {
              }
            }        
            catch (NoRuleException e) 
            {
              JOptionPane.showMessageDialog(topSplitPane_, e.getMessage()); 
              try 
              {
                writeStringInCurrentFile(fileBeforeSave);
              }
              catch (IOException e1)
              {             
              }
            }
          
            
          try 
          {          
            if(readFile(RulesContainer.currentRulesDir_+"/"+lipidClassName_ + "_" + lipidAdduct_ + ".frag.txt").contains("EMPTY"))            
              deleteFile(RulesContainer.currentRulesDir_+"/"+lipidClassName_ + "_" + lipidAdduct_ + ".frag.txt");
          }
          catch (IOException e) 
          {     
          }
        }      
      }
    });
    
    
    bottomSectionButtons_.add(showTotalDecisionButton);
    bottomSectionButtons_.add(saveRulesButton);
  }
  
  /**
   * Refreshes the container to get the new rules back - takes the default directory
   * @throws RulesException
   * @throws IOException
   * @throws SpectrummillParserException
   */
  private void setRuleClassIdentifier() throws RulesException, IOException, SpectrummillParserException 
  {
    ruleClassIdentifier_ = StaticUtils.getRuleName(lipidClassName_, lipidAdduct_); 
  }
  
  
  /**
   * Reading the current file
   * @param dateiName
   * @throws IOException
   */
  private String readFile(String dateiName) throws IOException
  {
    byte zeichen;    
    String text = "";    
    FileInputStream leseStrom = new FileInputStream(dateiName);
    do
    {
      zeichen = (byte)leseStrom.read();      
      text += (char)zeichen;
    } 
    while (zeichen !=-1);
    leseStrom.close();
    return text;
 }

  /**
   * Deletes a certain file
   * @param the path to the file
   */
  private void deleteFile(String path)
  {    
      File datei = new File(path);
      if (datei.exists()) 
      {
        datei.delete();        
      }
  }
  
    
  /**
   * Writes the buffer back in the file
   * @param input
   * @throws IOException
   */
  private void writeStringInCurrentFile(String input) throws IOException
  {
    String filename = lipidClassName_ + "_" + lipidAdduct_ + ".frag.txt";
    File file = new File("fragRules/"+filename);         
    if (!file.exists()) 
    {
      file.createNewFile();
    } 
    FileWriter fw = new FileWriter(file.getAbsoluteFile());
    BufferedWriter bw = new BufferedWriter(fw);
    if(input != null || input != "0")
      bw.write(input);
    bw.close();
  }
 
  
  /**
   * fetches the rules - first it tries to find the rules in the CACHE_DIR, and afterwards in the default location
   * @throws NoRuleException thrown if the rules are not there
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws IOException exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  private void getRulesFromCacheFirst() throws NoRuleException, RulesException, IOException, SpectrummillParserException
  {
    try{
      getRulesFromCertainDirectory(CACHE_DIR);
    } catch (NoRuleException nrx){
      try{
        getRules();
      } catch (NoRuleException nrx2){
        int reply = JOptionPane.showConfirmDialog(topSplitPane_, "This class does not exist. Create new one?", "Create new class" , JOptionPane.YES_NO_OPTION);
        if (reply == JOptionPane.YES_OPTION) 
        { 
          setEverythingNull();
          showFragmentTabs_ = false;
        }         
        else
          throw nrx2;
      }
    }
  }
  
  /**
   * Saves the rules from rules container in the global variables - takes the default directory
   * @throws NoRuleException catched by the settings button
   * @throws SpectrummillParserException 
   * @throws IOException 
   * @throws RulesException 
   */
  void getRules() throws NoRuleException, RulesException, IOException, SpectrummillParserException
  {      
    Double rtMaxDev = null;
    if (RulesContainer.getRetentionTimeMaxDeviation(ruleClassIdentifier_)!=null) rtMaxDev = new Double(RulesContainer.getRetentionTimeMaxDeviation(ruleClassIdentifier_));
    int faHydroxyRangeStart = -1;
    int faHydroxyRangeStop = -1;
    int lcbHydroxyRangeStart = -1;
    int lcbHydroxyRangeStop = -1;
    if (RulesContainer.getFaHydroxyRange(ruleClassIdentifier_)!=null && RulesContainer.getFaHydroxyRange(ruleClassIdentifier_).getStop()>0) {
      faHydroxyRangeStart = RulesContainer.getFaHydroxyRange(ruleClassIdentifier_).getStart();
      faHydroxyRangeStop = RulesContainer.getFaHydroxyRange(ruleClassIdentifier_).getStop();
    }
    if (RulesContainer.getLcbHydroxyRange(ruleClassIdentifier_)!=null && RulesContainer.getLcbHydroxyRange(ruleClassIdentifier_).getStop()>0) {
      lcbHydroxyRangeStart = RulesContainer.getLcbHydroxyRange(ruleClassIdentifier_).getStart();
      lcbHydroxyRangeStop = RulesContainer.getLcbHydroxyRange(ruleClassIdentifier_).getStop();
    }
    
    this.generalSettingsVO_ = new GeneralSettingsVO(new Integer(RulesContainer.getAmountOfChains(ruleClassIdentifier_)),
        new Integer(RulesContainer.getAmountOfAlkylChains(ruleClassIdentifier_)), new Integer(RulesContainer.getAmountOfAlkenylChains(ruleClassIdentifier_)),
        new Short(RulesContainer.getAmountOfLCBs(ruleClassIdentifier_)), RulesContainer.getAddChainPositions(ruleClassIdentifier_),
        RulesContainer.getChainlibrary(ruleClassIdentifier_), RulesContainer.getLcbLibrary(ruleClassIdentifier_), faHydroxyRangeStart,
        faHydroxyRangeStop, lcbHydroxyRangeStart, lcbHydroxyRangeStop, RulesContainer.getCAtomsFromNamePattern(ruleClassIdentifier_),
        RulesContainer.getDoubleBondsFromNamePattern(ruleClassIdentifier_), RulesContainer.isSingleChainIdentification(ruleClassIdentifier_),
        RulesContainer.getChainCutoffAsString(ruleClassIdentifier_), RulesContainer.getBasePeakCutoffAsString(ruleClassIdentifier_),
        RulesContainer.getSpectrumCoverageMinAsString(ruleClassIdentifier_), RulesContainer.isRtPostprocessing(ruleClassIdentifier_),
        RulesContainer.correctRtForParallelModel(ruleClassIdentifier_), rtMaxDev,
        RulesContainer.getMSIdentificationOrder(ruleClassIdentifier_));
    
      this.headFragmentRules_ = RulesContainer.getHeadFragmentRules(ruleClassIdentifier_);
      this.headIntensityRules_ = RulesContainer.getHeadIntensityRules(ruleClassIdentifier_);
      this.chainFragmentRules_ = RulesContainer.getChainFragmentRules(ruleClassIdentifier_);
      this.chainIntensityRules_ = RulesContainer.getChainIntensityRules(ruleClassIdentifier_);      
      this.positionIntensityRules_  = RulesContainer.getPositionIntensityRules(ruleClassIdentifier_);    
  }
  
/**
 * Saves the rules from rules container in the global variables from a certain directory
 * @param fileDir the caach directory for the rule txt file
 * @throws NoRuleException
 * @throws RulesException
 * @throws IOException
 * @throws SpectrummillParserException
 */
  void getRulesFromCertainDirectory(String fileDir) throws NoRuleException, RulesException, IOException, SpectrummillParserException
  {         
    Double rtMaxDev = null;
    if (RulesContainer.getRetentionTimeMaxDeviation(ruleClassIdentifier_,fileDir)!=null) rtMaxDev = new Double(RulesContainer.getRetentionTimeMaxDeviation(ruleClassIdentifier_,fileDir));
    int faHydroxyRangeStart = -1;
    int faHydroxyRangeStop = -1;
    int lcbHydroxyRangeStart = -1;
    int lcbHydroxyRangeStop = -1;
    if (RulesContainer.getFaHydroxyRange(ruleClassIdentifier_,fileDir)!=null && RulesContainer.getFaHydroxyRange(ruleClassIdentifier_,fileDir).getStop()>0) {
      faHydroxyRangeStart = RulesContainer.getFaHydroxyRange(ruleClassIdentifier_,fileDir).getStart();
      faHydroxyRangeStop = RulesContainer.getFaHydroxyRange(ruleClassIdentifier_,fileDir).getStop();
    }
    if (RulesContainer.getLcbHydroxyRange(ruleClassIdentifier_,fileDir)!=null && RulesContainer.getLcbHydroxyRange(ruleClassIdentifier_,fileDir).getStop()>0) {
      lcbHydroxyRangeStart = RulesContainer.getLcbHydroxyRange(ruleClassIdentifier_,fileDir).getStart();
      lcbHydroxyRangeStop = RulesContainer.getLcbHydroxyRange(ruleClassIdentifier_,fileDir).getStop();
    }

    this.generalSettingsVO_ = new GeneralSettingsVO(new Integer(RulesContainer.getAmountOfChains(ruleClassIdentifier_,fileDir)),
        new Integer(RulesContainer.getAmountOfAlkylChains(ruleClassIdentifier_,fileDir)), new Integer(RulesContainer.getAmountOfAlkenylChains(ruleClassIdentifier_,fileDir)),
        new Short(RulesContainer.getAmountOfLCBs(ruleClassIdentifier_,fileDir)),RulesContainer.getAddChainPositions(ruleClassIdentifier_,fileDir),
        RulesContainer.getChainlibrary(ruleClassIdentifier_,fileDir), RulesContainer.getLcbLibrary(ruleClassIdentifier_,fileDir),
        faHydroxyRangeStart, faHydroxyRangeStop, lcbHydroxyRangeStart, lcbHydroxyRangeStop,RulesContainer.getCAtomsFromNamePattern(ruleClassIdentifier_,fileDir), 
        RulesContainer.getDoubleBondsFromNamePattern(ruleClassIdentifier_,fileDir), RulesContainer.isSingleChainIdentification(ruleClassIdentifier_,fileDir),
        RulesContainer.getChainCutoffAsString(ruleClassIdentifier_,fileDir), RulesContainer.getBasePeakCutoffAsString(ruleClassIdentifier_,fileDir),
        RulesContainer.getSpectrumCoverageMinAsString(ruleClassIdentifier_,fileDir), RulesContainer.isRtPostprocessing(ruleClassIdentifier_,fileDir),
        RulesContainer.correctRtForParallelModel(ruleClassIdentifier_,fileDir), rtMaxDev,
        RulesContainer.getMSIdentificationOrder(ruleClassIdentifier_,fileDir));

      this.headFragmentRules_ = RulesContainer.getHeadFragmentRules(ruleClassIdentifier_, fileDir);
      this.headIntensityRules_ = RulesContainer.getHeadIntensityRules(ruleClassIdentifier_, fileDir);
      this.chainFragmentRules_ = RulesContainer.getChainFragmentRules(ruleClassIdentifier_, fileDir);
      this.chainIntensityRules_ = RulesContainer.getChainIntensityRules(ruleClassIdentifier_, fileDir);      
      this.positionIntensityRules_  = RulesContainer.getPositionIntensityRules(ruleClassIdentifier_, fileDir);    
  }
  
  /**
   * Refreshes the middleTab section if a rule has changed but ignores the current settings in the general tab
   * @param tabIndex Index of the current tab
  */
  void refreshMiddleWithoutCurrentGenerals(int tabIndex)
  {
    topSplitPane_.remove( middleTabSection_ );
    topSplitPane_.validate();
    ruleTabsSection();
    middleTabSection_.setSelectedIndex(tabIndex);
    topSplitPane_.setBottomComponent(middleTabSection_);  
  }

  /**
   * Refreshes the middleTab section if a rule has changed
   * @param tabIndex Index of the current tab
   */
  void refreshMiddle(int tabIndex)
  {
    this.generalSettingsVO_ = this.generalSettings_.getValues();
    
    topSplitPane_.remove( middleTabSection_ );
    topSplitPane_.validate();
    ruleTabsSection();
    middleTabSection_.setSelectedIndex(tabIndex);
    topSplitPane_.setBottomComponent(middleTabSection_);  
  }
  
  
  /**
   * saves the rules into a directory and updates the RulesContainer
   * @param dir the directory the rules shall be saved into
   * @throws IOException exception if there is something wrong about the file
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   * @throws NoRuleException thrown if the rules are not there
   */
  public void saveRules(String dir) throws IOException, RulesException, SpectrummillParserException, NoRuleException{
    writeRules(dir);
    if(dir == "default")
    {
      RulesContainer.clearCache();
      getRulesFromCertainDirectory(RulesContainer.currentRulesDir_);   
//      paintNewSpectra(false);
    }else{
      RulesContainer.clearCache(dir);
    }
    deletePossibleDummy();  

  }
  
  /**
   * Prints all the rules in the special class text file
   * @throws IOException exception if there is something wrong about the file
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   * @throws NoRuleException thrown if the rules are not there
   */
  public void writeRules(String dir) throws IOException, RulesException, SpectrummillParserException, NoRuleException
  {  
    FragRuleParser.writeRules(dir, lipidClassName_, lipidAdduct_, generalSettingsVO_, headFragments_, this.headIntensityRules_,
        sortedChainFragments_, chainIntensityRules_, positionIntensityRules_);
  }
  
  
  
  public void saveDummy(String dir) throws IOException, RulesException, SpectrummillParserException, NoRuleException
  {  
    File file;
    String filename = lipidClassName_ + "_" + lipidAdduct_ + ".frag.txt";    
    if(dir == "default")
    {
      file = new File("fragRules/"+filename);   
    }
    else
    {
      file = new File(dir + "/"+filename); 
    }
    if (!file.exists())      
    {
      file.createNewFile();
    } 
    FileWriter fw = new FileWriter(file.getAbsoluteFile());
    BufferedWriter bw = new BufferedWriter(fw);

    bw.write("[GENERAL]\n");        
    bw.write("AmountOfChains="+ "0" + "\n");             
    bw.write("ChainLibrary="+ "fattyAcidChains.xlsx" + "\n");        
    bw.write("CAtomsFromName="+ "(ENTER)" + "\n");        
    bw.write("DoubleBondsFromName="+ "(ENTER)" + "\n");       
    
    bw.write("BasePeakCutoff="+ "0.0%" + "\n");   
    bw.write("SpectrumCoverage="+ "0.0%" + "\n");    
        
    bw.write("RetentionTimePostprocessing=" + "false" + "\n");    
    bw.write("RetentionTimeParallelSeries=" + "false" + "\n");
    bw.write("\n");
    bw.write("[HEAD]\n");
    bw.write("!FRAGMENTS\n");        

    FragmentRuleVO ruleVO = new FragmentRuleVO("EMPTY","CCCC",1, 2, (short)0, headFragmentRules_, chainFragmentRules_, elementParser_);
    headFragmentRules_.put("EMPTY", ruleVO);
    headFragments_ = getSortedVector(headFragmentRules_, TYPE_HEAD);  
    String mandatory;    
    if(headFragments_.size() != 0)
    {
    for(int i=0; i<headFragments_.size(); i++ )
    {
      mandatory = "false";                
      if(headFragments_.elementAt(i).isMandatory()==FragmentRuleVO.MANDATORY_TRUE)
      {
        mandatory = "true";
      }          
      bw.write("Name=" + headFragments_.elementAt(i).getName() + "\t" + "Formula=" + 
        headFragments_.elementAt(i).getFormula() + "\t" + "Charge=" + 
      Integer.toString(headFragments_.elementAt(i).getCharge()) + "\t"+ "MSLevel=" + 
      Integer.toString(headFragments_.elementAt(i).getMsLevel()) + "\t" + "mandatory=" + mandatory + "\n");     
    }  
    bw.close();
    }
  }
  
    
  /**
   * Return the name string mandatory
   * @param mandatoryInput
   * @return mandatory as a string
   */
  String getStringMandatory(boolean mandatoryInput)
  {  
    String mandatory = "false";           
    if(mandatoryInput)
    {
      mandatory = "true";
    }   
    return mandatory;
  }
  
  /**
   * Sorts the FragmentVOs in a hashtable and puts it into a vector.
   * @param input
   * @return the sorted vector
   */
  private Vector<FragmentRuleVO> getSortedVector(Hashtable<String,FragmentRuleVO> input, int type)
  {     
    Hashtable<String,FragmentRuleVO> cache = new Hashtable<String,FragmentRuleVO>(input);
    Vector<FragmentRuleVO> output = new Vector<FragmentRuleVO>();
    Hashtable<String,FragmentRuleVO> outputHash = new Hashtable<String,FragmentRuleVO>();
    int before = Integer.MIN_VALUE;
    while (cache.size()>0 && before!=cache.size()){
      Vector<String> usedKeys = new Vector<String>();
      for (String key : cache.keySet()){
        FragmentRuleVO rule = cache.get(key);
        try{
          if (type==TYPE_HEAD){
            new FragmentRuleVO(rule.getName(),rule.getFormula(),rule.getCharge(),rule.getMsLevel(),rule.isMandatory(),
                outputHash, new Hashtable<String,FragmentRuleVO>(), elementParser_);
          } else if (type==TYPE_CHAIN){
            new FragmentRuleVO(rule.getName(),rule.getFormula(),rule.getCharge(),rule.getMsLevel(),rule.isMandatory(),
                headFragmentRules_, outputHash, elementParser_);            
          }
          outputHash.put(key, rule);
          output.add(rule);
          usedKeys.add(key);
        // if there is an undefined Fragment, the RulesException is thrown  
        } catch (RulesException rx){}
      }
      before = cache.size();
      //concurrent access on the hashtable is not allowed in a for loop, thus
      // I have to remove the found ones separately
      for (String key : usedKeys){
        cache.remove(key);
      }
    }
    for (String key : cache.keySet()){
      output.add(cache.get(key));
    }
    return output;
  }  

  
  /**
   * Returns the current combo setting
   */
  private short getFragMandatory(int inputNumber)
  {
    if(inputNumber == 0)
    {
      return 1;          
    }
    else
    {   
      return 0;
    }
  }

  /**
   * Returns the current combo setting
   */
  private boolean getIntensityMandatory(int inputNumber)
  {
    if(inputNumber == 0)
    {
      return true;          
    }
    else
    {   
      return false;
    }
  }

  /**
   * Builds the msn analyzer and returns the current debug info.
   * @return Current debug info
   * @throws CgException 
   * @throws SpectrummillParserException 
   * @throws IOException 
   * @throws RulesException 
   * @throws HydroxylationEncodingException thrown if hydroxylation the encoding does not exist
   * @throws ChemicalFormulaException thrown if there is something wrong with the formula
   * @throws LipidCombinameEncodingException thrown when a lipid combi id (containing type and OH number) cannot be decoded
   */
  MSnDebugVO getCurrentDebugInfo() throws RulesException, IOException, SpectrummillParserException, CgException, HydroxylationEncodingException, ChemicalFormulaException, LipidCombinameEncodingException
  {    
    MSnDebugVO debugInfo = null;  
    try 
    {
      msnAnalyzer_ = updateMSnAnalyzerToCurrentSettings();
    }
    catch (NoRuleException e) 
    {
      e.printStackTrace();
    }
    debugInfo = msnAnalyzer_.getDebugInfo(); 
    return debugInfo;
  }
  
  
  /**
   * Returns a string with the labeled probes of each fragment
   * @param detectedFragments
   * @return
   * @throws CgException 
   * @throws SpectrummillParserException 
   * @throws IOException 
   * @throws RulesException 
   * @throws HydroxylationEncodingException thrown if hydroxylation the encoding does not exist
   * @throws ChemicalFormulaException thrown if there is something wrong with the formula
   * @throws LipidCombinameEncodingException thrown when a lipid combi id (containing type and OH number) cannot be decoded
   */
  private String showFragmentTabDetailsString(Hashtable<String,CgProbe> detectedFragments) throws RulesException, IOException, SpectrummillParserException, CgException, HydroxylationEncodingException, ChemicalFormulaException, LipidCombinameEncodingException
  {      
      MSnDebugVO debugInfo = getCurrentDebugInfo();      
    
      String labeling = "Area";              
      
      Enumeration<String> headFragmentEnum = detectedFragments.keys();
      
      String labelInput = "<br>Fragments:<br>";      
      while (headFragmentEnum.hasMoreElements()) 
      {
        String currentKey = (String)headFragmentEnum.nextElement();   
        CgProbe inhalt = (CgProbe) detectedFragments.get(currentKey); 
        String[] fragmentRuleResults = inhalt.toString().split(";");
        String[] labels = labeling.toString().split(";");
        
        labelInput = labelInput + "<br><font color='red'>" + currentKey + "</font><br>";
        if(fragmentRuleResults.length != 0)
        {      
          labelInput = labelInput + labels[0] +  " = " + fragmentRuleResults[0] + "<br>";         
        }
      }         
               
          
      Hashtable<String, Integer> discHeads = debugInfo.getDiscardedHeadGroupFragments();
      labelInput = labelInput + "<br>Not found fragments: " + discHeads.size();
      for (String fragment: discHeads.keySet())
      {
        labelInput = labelInput + "<br><br>Head-Discarded: "+fragment+"\t"+ "because " + printStatus(discHeads.get(fragment));
      }
      
      return labelInput;    
  }
  
  /**
   * Returns a string with the used head intensity rules
   * @param fulfilledIntensityRules
   * @return
   * @throws CgException 
   * @throws SpectrummillParserException 
   * @throws IOException 
   * @throws RulesException 
   * @throws HydroxylationEncodingException thrown if hydroxylation the encoding does not exist
   * @throws ChemicalFormulaException thrown if there is something wrong with the formula
   * @throws LipidCombinameEncodingException thrown when a lipid combi id (containing type and OH number) cannot be decoded
   */
  private String showHeadIntensityTabDetailsString(Hashtable<String,IntensityRuleVO> fulfilledIntensityRules) throws RulesException, IOException, SpectrummillParserException, CgException, HydroxylationEncodingException, ChemicalFormulaException, LipidCombinameEncodingException
  {   
      MSnDebugVO debugInfo = getCurrentDebugInfo();   
      String labelInput = "Intensity rules that were fulfilled by the head group fragments:<br><br>";      
      
      Enumeration<String> headIntensityEnum = fulfilledIntensityRules.keys(); 
      
      int counter = 0;
      while (headIntensityEnum.hasMoreElements()) 
      {
        counter++;
        String currentKey = (String)headIntensityEnum.nextElement();   
        IntensityRuleVO inhalt = (IntensityRuleVO) fulfilledIntensityRules.get(currentKey); 
        
        labelInput = labelInput + counter + ". " + "Equation= " + inhalt.getRuleIdentifier() + " " + "Mandatory= " + inhalt.isMandatory() + "<br>";                
      }    
      
      Hashtable<String,IntensityRuleVO> violHeadRules = debugInfo.getViolatedHeadRules();
      labelInput = labelInput + "<br><br>Violated head intensity rules: " + violHeadRules.size();  
      
      return labelInput;   
  }
  
  /**
   * Returns a string of the intensity rules used
   * @param fulfilledIntensityRules
   * @return
   * @throws LipidCombinameEncodingException thrown when a lipid combi id (containing type and OH number) cannot be decoded
   */
  private String showChainIntensityTabDetailsString(Hashtable<String,IntensityChainVO> chainIntensities) throws LipidCombinameEncodingException
  {  
    int counter = 0;
    String returnString = "";
    Enumeration<String> chainIntensitiesEnum = chainIntensities.keys();       
    while (chainIntensitiesEnum.hasMoreElements()) 
    {  
      counter++;
      String currentKey = (String)chainIntensitiesEnum.nextElement();  
      IntensityChainVO intensityVo = (IntensityChainVO) chainIntensities.get(currentKey);      
      returnString = returnString + counter + ". " + intensityVo.getReadableRuleInterpretation(Settings.getFaHydroxyEncoding(),Settings.getLcbHydroxyEncoding());
    }
    return returnString;
  }
  
  /**
   * Prints out a dialog of the possible positions of the chains
   * @param result the obtained result
   * @param positionRecommendations1
   * @throws CgException 
   * @throws SpectrummillParserException 
   * @throws IOException 
   * @throws RulesException 
   * @throws HydroxylationEncodingException thrown if hydroxylation the encoding does not exist
   * @throws ChemicalFormulaException thrown if there is something wrong with the formula
   * @throws LipidCombinameEncodingException thrown when a lipid combi id (containing type and OH number) cannot be decoded
   */
  public void printPositionRecomandation (LipidParameterSet result, Hashtable<String,Hashtable<Integer,Integer>> positionRecommendations1) throws RulesException, IOException, SpectrummillParserException, CgException, HydroxylationEncodingException, ChemicalFormulaException, LipidCombinameEncodingException
  {   
    deleteDetailsBoxes();   
    
    MSnDebugVO debugInfo = getCurrentDebugInfo();
    
    String toWrite ="";    
    toWrite = toWrite + "<br>Possible Positions:<br>" + getResultWriteString(result,false);
        
    Hashtable<String,Hashtable<String,IntensityRuleVO>> unfulfilledPosRules = debugInfo.getUnfulfilledPositionRules();
    toWrite = toWrite + "<br><br>Unfulfilled rules: "+unfulfilledPosRules.size()+"<br/>";
    for (String combiName : unfulfilledPosRules.keySet())
    {
    	String humanReadable = StaticUtils.getHumanReadableCombiName(combiName,Settings.getFaHydroxyEncoding(),Settings.getLcbHydroxyEncoding());
      Hashtable<String,IntensityRuleVO> unfulfilled = unfulfilledPosRules.get(combiName);
      for (IntensityRuleVO ruleVO : unfulfilled.values())
      {
        toWrite = toWrite + humanReadable+":\tNOT "+ruleVO.getReadableRuleInterpretation(Settings.getFaHydroxyEncoding(),Settings.getLcbHydroxyEncoding())+"<br/>";
      }
    }
    toWrite += "<br/>";
    Hashtable<String,Vector<Vector<IntensityRuleVO>>> contradictingPositionRules  = debugInfo.getContradictingPositionRules();
    toWrite = toWrite + "<br>Any contradicting position rules: "+contradictingPositionRules.size() + "<br>";
    for (String combiName : contradictingPositionRules.keySet())
    {
      Vector<Vector<IntensityRuleVO>> rulePairs = contradictingPositionRules.get(combiName);
      for (Vector<IntensityRuleVO> rulePair : rulePairs)
      {
        IntensityRuleVO rule1 = rulePair.get(0);
        IntensityRuleVO rule2 = rulePair.get(1);
        toWrite = toWrite + combiName+":\t"+rule1.getReadableRuleInterpretation(Settings.getFaHydroxyEncoding(),Settings.getLcbHydroxyEncoding())+"\t!=\t"
        +rule2.getReadableRuleInterpretation(Settings.getFaHydroxyEncoding(),Settings.getLcbHydroxyEncoding());
      }
    }
    
    
    JLabel lbl = new JLabel("<html>"+toWrite+"<br></html>");
    lbl.setVerticalAlignment(JLabel.TOP);  
    
    showHeadDetails_ = new JDialog();
    showHeadDetails_.setTitle("Show Details");
    showHeadDetails_.setSize(600,300);
    
    JScrollPane headRule = new JScrollPane(lbl);       
    showHeadDetails_.add(headRule);
    showHeadDetails_.setVisible(true); 
  }
  
  /**
   * Delets the details box
   */
  public void deleteDetailsBoxes()
  {
    if(showHeadDetails_ != null)
    {
      showHeadDetails_.dispose();
    }
  }
  
  /**
   * Gets the Msn Analyzer hashes of the rules used and prints it together
   * @param detectedFragments
   * @param fulfilledIntensityRules
   * @param chainFragments
   * @param chainIntensities
   * @throws CgException 
   * @throws SpectrummillParserException 
   * @throws IOException 
   * @throws RulesException 
   * @throws HydroxylationEncodingException thrown if hydroxylation the encoding does not exist
   * @throws ChemicalFormulaException thrown if there is something wrong with the formula
   * @throws LipidCombinameEncodingException thrown when a lipid combi id (containing type and OH number) cannot be decoded
   */
  private void printInDialogShowDetailsBox(Hashtable<String,CgProbe> detectedFragments, 
      Hashtable<String,IntensityRuleVO> fulfilledIntensityRules, 
      Hashtable<String,Hashtable<String,CgProbe>> chainFragments, Hashtable<String,Hashtable<String,IntensityChainVO>> chainIntensities) throws RulesException, IOException, SpectrummillParserException, CgException, HydroxylationEncodingException, ChemicalFormulaException, LipidCombinameEncodingException
  {
    try{
    deleteDetailsBoxes();
    
    MSnDebugVO debugInfo = getCurrentDebugInfo();   
    String toWrite = "";
    if(detectedFragments != null)
    {
      toWrite = showFragmentTabDetailsString(detectedFragments) + 
      "<br><br>" + showHeadIntensityTabDetailsString(fulfilledIntensityRules);
    }
    else
    {
      ArrayList<String> chainNamesInternalRepresentation = new ArrayList<String>(chainFragments.keySet());
      toWrite = toWrite + "<br>Detected Chains<br><br>";
      for (String name : chainNamesInternalRepresentation)
      {
      	String humanReadable = StaticUtils.getHumanReadableCombiName(name,Settings.getFaHydroxyEncoding(),Settings.getLcbHydroxyEncoding());
      	toWrite = toWrite + "Chain " + humanReadable + "<br>"; 
      }
      toWrite = toWrite + "<br><br>";  
      
      for (String name : chainNamesInternalRepresentation)
      {
      	String humanReadable = StaticUtils.getHumanReadableCombiName(name,Settings.getFaHydroxyEncoding(),Settings.getLcbHydroxyEncoding());
      	toWrite = toWrite + "Details Chain " + humanReadable + "<br>";
        Hashtable<String,CgProbe> chainFragment = (Hashtable<String,CgProbe>) chainFragments.get(name);
        toWrite = toWrite + showFragmentTabDetailsString(chainFragment) + "<br><br>"; 
      }
            
      toWrite = toWrite + "<br><br>Intensity rules that were fulfilled by the chain fragments:<br><br>";
      Enumeration<String> chainIntensityEnum = chainIntensities.keys();         
      while (chainIntensityEnum.hasMoreElements()) 
      { 
      	String currentKey = (String)chainIntensityEnum.nextElement();
        Hashtable<String,IntensityChainVO> fulfilledChainIntensityRules = (Hashtable<String,IntensityChainVO>) chainIntensities.get(currentKey);       
        toWrite = toWrite + showChainIntensityTabDetailsString(fulfilledChainIntensityRules) + "<br>";
      }
      toWrite = toWrite + "<br><br><br>";
      
      Hashtable<String,Hashtable<String,Object>> violChainRules = debugInfo.getViolatedChainRules();
      toWrite = toWrite + "Violated chain intensity rules: " + violChainRules.size() + "<br>";
      for (String faName : violChainRules.keySet())
      {
      	String humanReadable = StaticUtils.getHumanReadableCombiName(faName,Settings.getFaHydroxyEncoding(),Settings.getLcbHydroxyEncoding());
        Hashtable<String,Object> violRule = violChainRules.get(faName);
        for (String ruleName : violRule.keySet())
        {
          Object rule = violRule.get(ruleName);
          if (rule instanceof Integer)
          {
            int status = (Integer)rule;
            toWrite = toWrite + "<br>Chain Discarded: " + humanReadable+"\t"+ruleName+"\t"+ "because "+ printStatus(status);
          }
          else if (rule instanceof IntensityRuleVO)
          {
	          IntensityRuleVO ruleVO = (IntensityRuleVO)rule;
	          toWrite = toWrite + "<br>Violated chain rule: "+ humanReadable+"\t"+ruleVO.getReadableRuleInterpretation(Settings.getFaHydroxyEncoding(),Settings.getLcbHydroxyEncoding());
          }
        }
      }
      Hashtable<String,Integer> violCombis = debugInfo.getViolatedCombinations();
      toWrite = toWrite + "<br><br>Violated combinations: "+violCombis.size() + "<br>";
      for (String combi : violCombis.keySet())
      {
      	String humanReadable = StaticUtils.getHumanReadableCombiName(combi,Settings.getFaHydroxyEncoding(),Settings.getLcbHydroxyEncoding());
        toWrite = toWrite + "<br>Violated combi: "+humanReadable+"\t"+violCombis.get(combi);
      }
      toWrite = toWrite + "<br>Spectrum sufficiently covered: "+debugInfo.isSpectrumCoverageFulfilled();
    }  
    
      
    JLabel lbl = new JLabel("<html>"+toWrite+"<br></html>");
    lbl.setVerticalAlignment(JLabel.TOP);  
    
    showHeadDetails_ = new JDialog();
    showHeadDetails_.setTitle("Show Details");
    showHeadDetails_.setSize(600,800);
    
    JScrollPane headRule = new JScrollPane(lbl);       
    showHeadDetails_.add(headRule);
    showHeadDetails_.setVisible(true); 
    } catch (RulesException rx){
      if (rx.getMessage().endsWith(FragRuleParser.NO_HEAD_AND_CHAINS_SECTION)){
        new WarningMessage(new JFrame(), "Error", "There are no rules entered!");
      } else throw rx;
    }
  }   
  
  /**
   * Prints the status of the debugVO
   * @param status
   * @return the string status
   */
  String printStatus(int status)
  {
    String returnStatus = null;
    
    if(status == 0)
      returnStatus = "the reason for the discard is unknown"; 
    if(status == 1)
      return returnStatus = "there is no peak at this m/z value";    
    if(status == 2)
      return returnStatus = "the peak is lower than the base peak intensity cutoff";    
    if(status == 3)
      return returnStatus = "the found chain does not have any matching partner chains to match the total number of carbon atoms and double bonds";       
    if(status == 4)
      return returnStatus = "the found chain combination is of minor intensity - below the defined threshold";
    
    return returnStatus;
  }
  
  
  public void deleteDependingRules(String name, int type){    
    Vector<String> names = this.deleteDependendFragment(name, type);
    this.deleteDependendEquation(names, type);
  }
  
  /**
   * Delets a dependend fragment
   * @param deletedName
   * @param type 1 for head, 2 for chain
   */
  public Vector<String> deleteDependendFragment(String deletedName, int type)
  {
    Vector<String> deleted = new Vector<String>();
    deleted.add(deletedName);
    Hashtable<String,FragmentRuleVO> allowedRules = new Hashtable<String,FragmentRuleVO>();
    if (type==TYPE_HEAD){
      Hashtable<String,FragmentRuleVO> allowedHeadRules = new Hashtable<String,FragmentRuleVO>();
      Vector<FragmentRuleVO> sortedFragments = this.getSortedVector(headFragmentRules_, TYPE_HEAD);
      for (FragmentRuleVO ruleVO : sortedFragments){
        if (ruleVO.containsOnlyAllowedFragments(allowedRules)){
          allowedHeadRules.put(ruleVO.getName(), ruleVO);
          allowedRules.put(ruleVO.getName(), ruleVO);
        }else{
          deleted.add(ruleVO.getName());
        }
      }
      headFragmentRules_ = allowedHeadRules;
    }else{
      allowedRules = new Hashtable<String,FragmentRuleVO>(headFragmentRules_);
    }
    if (type<TYPE_POSITION){
      Hashtable<String,FragmentRuleVO> allowedChainRules = new Hashtable<String,FragmentRuleVO>();
      Vector<FragmentRuleVO> sortedFragments = this.getSortedVector(chainFragmentRules_, TYPE_CHAIN);
      for (FragmentRuleVO ruleVO : sortedFragments){
        if (ruleVO.containsOnlyAllowedFragments(allowedRules)){
          allowedChainRules.put(ruleVO.getName(), ruleVO);
          allowedRules.put(ruleVO.getName(), ruleVO);
        }else{
          deleted.add(ruleVO.getName());
        }
      }
      chainFragmentRules_ = allowedChainRules;
    }
    return deleted;
  }
  

  /**
   * Looks for equations that depend on the fragments and deletes it
   * @param name as the current key of the hashtable
   * @param type 1...head, 2... chain
   */
  public void deleteDependendEquation(Vector<String> names, int type)
  {
    if (type==TYPE_HEAD) headIntensityRules_ = removeDependentIntensityRule(headIntensityRules_, names);
    if (type<TYPE_POSITION) chainIntensityRules_ = removeDependentIntensityRule(chainIntensityRules_, names);
    positionIntensityRules_ = removeDependentIntensityRule(positionIntensityRules_, names);
  }
  
  private Vector<IntensityRuleVO> removeDependentIntensityRule(Vector<IntensityRuleVO> rules, Vector<String> fragmentNames){
    Vector<IntensityRuleVO> passingRules = new Vector<IntensityRuleVO>();
    for (IntensityRuleVO rule : rules){
      boolean hasWrongFragment = false;
      for (String name : fragmentNames){
        if (rule.containsFragment(name)) hasWrongFragment = true;
      }
      if (!hasWrongFragment) passingRules.add(rule);
    }
    return passingRules;
  }
  
  
  /**
   * Prints a pop up box with the total decision result
   */
  public void printTotalDecisionInBox()
  {
    try{
    deleteDetailsBoxes();
    
    msnAnalyzer_ = null;    
    try 
    {
      msnAnalyzer_ = updateMSnAnalyzerToCurrentSettings();
    }
    catch (IOException e1) 
    {
      e1.printStackTrace();
    }
    catch (SpectrummillParserException e1) 
    {
      e1.printStackTrace();
    }
    catch (CgException e1) 
    {
      e1.printStackTrace();
    }     
    catch (NoRuleException e) 
    {
    }
    catch (HydroxylationEncodingException | ChemicalFormulaException | LipidCombinameEncodingException e) {
      e.printStackTrace();
    }  
    LipidParameterSet result = msnAnalyzer_.getResult();  
    
    String toWrite = "";
    toWrite = toWrite + "<br><font color='red'>Result:<br>" + getResultWriteString(result,true) + "</font>";    
    toWrite = toWrite + "<br><br>Double Bonds: " + result.getDoubleBonds();
    toWrite = toWrite + "<br><br>Modification Name: " + result.getModificationName();
    toWrite = toWrite + "<br><br>Analyte Formula: " + result.getAnalyteFormula();
    toWrite = toWrite + "<br><br>Modification Formula: " + result.getModificationFormula();
    toWrite = toWrite + "<br><br>Chemical Formula: " + result.getChemicalFormula();
    toWrite = toWrite + "<br><br>Charge: " + result.getCharge();
    toWrite = toWrite + "<br><br>Rt: " + result.getRt();    
    
    
    JLabel lbl = new JLabel("<html>"+toWrite+"<br></html>");
    lbl.setVerticalAlignment(JLabel.TOP);  
    
    showHeadDetails_ = new JDialog();
    showHeadDetails_.setTitle("Total Decision");
    showHeadDetails_.setSize(600,800);
    
    JScrollPane headRule = new JScrollPane(lbl);       
    showHeadDetails_.add(headRule);
    showHeadDetails_.setVisible(true); 
    } catch (RulesException rx){
      if (rx.getMessage().endsWith(FragRuleParser.NO_HEAD_AND_CHAINS_SECTION) || !this.showFragmentTabs_){
        new WarningMessage(new JFrame(), "Error", "There are no rules entered!");
      } else rx.printStackTrace();;
    }  
    }  
  
  /**
   * Checks if a head fragment is present in the spectra
   * @throws IOException
   * @throws RulesException
   * @throws SpectrummillParserException
   * @throws NoRuleException
   * @throws CgException
   * @throws HydroxylationEncodingException thrown if hydroxylation the encoding does not exist
   * @throws ChemicalFormulaException thrown if there is something wrong with the formula
   * @throws LipidCombinameEncodingException thrown when a lipid combi id (containing type and OH number) cannot be decoded
   */
  public void checkToAddHeadFragment() throws IOException, RulesException, SpectrummillParserException, NoRuleException, CgException, HydroxylationEncodingException, ChemicalFormulaException, LipidCombinameEncodingException 
  { 
    //if (!this.checkGeneralEntries(true)) return;
    int before = 0;
    int after = 0;
    
    headRuleTab_.remove( fragRuleOK_ );   
    fragRuleOK_ = new JLabel( new ImageIcon( redCrossPicture_ ));
    gridBayLayout_.setConstraints(fragRuleOK_, new GridBagConstraints (0, headFragments_.size()+1, 1, 1, 0, 0, 
    GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0 ));
    headRuleTab_.add( fragRuleOK_ ); 
    headRuleTab_.repaint();
    headRuleTab_.validate();
    //The container has to be set back
    paintNewSpectra(false); 
    before = headFragmentCounter(); 
    FragmentRuleVO ruleVO = null;
    if(headFragmentRules_.containsKey(fragRuleNameField_.getText())) return;
    if(checkHeadFragment())
	    {       
	      ruleVO = new FragmentRuleVO(fragRuleNameField_.getText(),fragRuleFormulaField_.getText(),Integer.parseInt(fragRuleChargeField_.getText()), 
	      Integer.parseInt(fragRuleMsLevelField_.getText()), getFragMandatory(fragRuleMandatoryCombo_.getSelectedIndex()), headFragmentRules_, chainFragmentRules_, elementParser_);        
	      headFragmentRules_.put(fragRuleNameField_.getText(), ruleVO); 
	      headFragments_ = getSortedVector(headFragmentRules_, TYPE_HEAD);      
	    } else return; 
	    paintNewSpectra(false);      
	    after = headFragmentCounter();             
	    
	    if(after ==  before+1)
	    {      
	      headRuleTab_.remove( fragRuleOK_ );   
	      fragRuleOK_ = new JLabel( new ImageIcon( resource1_ ));
	      gridBayLayout_.setConstraints(fragRuleOK_, new GridBagConstraints (0, headFragments_.size(), 1, 1, 0, 0, 
	      GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0 ));
	      headRuleTab_.add( fragRuleOK_ ); 
	      headRuleTab_.repaint();
	      headRuleTab_.validate(); 
	    }    
	    
	    if(checkHeadFragment())
	    {  
	      headFragmentRules_.remove(fragRuleNameField_.getText());       
	      headFragments_ = getSortedVector(headFragmentRules_, TYPE_HEAD); 
	    }
  }
  
  /**
   * Checks if there is any problem while building the head fragmentVO
   * @return true for no problems, false for an error happened
   */
  public boolean checkHeadFragment()
  {
    boolean check = true;
    try 
    {
      new FragmentRuleVO(fragRuleNameField_.getText(),fragRuleFormulaField_.getText(),Integer.parseInt(fragRuleChargeField_.getText()), 
      Integer.parseInt(fragRuleMsLevelField_.getText()), getFragMandatory(fragRuleMandatoryCombo_.getSelectedIndex()), headFragmentRules_, chainFragmentRules_, elementParser_);
    }
    catch (NumberFormatException e) 
    {
      check = false;
    }
    catch (RulesException e) 
    {
      check = false;
    }
    if (fragRuleNameField_.getText()==null || fragRuleNameField_.getText().length()==0 ||
        fragRuleFormulaField_.getText()==null || fragRuleFormulaField_.getText().length()==0 ||
        fragRuleChargeField_.getText()==null || fragRuleChargeField_.getText().length()==0 ||
        fragRuleMsLevelField_.getText()==null || fragRuleMsLevelField_.getText().length()==0)
      check = false;
    if (!FragmentVerifier.checkFragmentNameInput(fragRuleNameField_)) check = false;
    if (!FragmentVerifier.checkFragmentFormulaInput(fragRuleFormulaField_)) check = false;
    return check;
  }
  
  /**
   * Counts the detected head fragments
   * @return the number of detected head fragments
   * @throws RulesException
   * @throws IOException
   * @throws SpectrummillParserException
   * @throws CgException
   * @throws HydroxylationEncodingException thrown if hydroxylation the encoding does not exist
   * @throws ChemicalFormulaException thrown if there is something wrong with the formula
   * @throws LipidCombinameEncodingException thrown when a lipid combi id (containing type and OH number) cannot be decoded
   */
  public int headFragmentCounter() throws RulesException, IOException, SpectrummillParserException, CgException, HydroxylationEncodingException, ChemicalFormulaException, LipidCombinameEncodingException
  {
    int numberOfDetectedFragments = 0;    
    try 
    {
      msnAnalyzer_ = updateMSnAnalyzerToCurrentSettings();
      Hashtable<String,CgProbe> detectedFragments = msnAnalyzer_.getHeadGroupFragments(); 
      numberOfDetectedFragments = detectedFragments.size();
      return numberOfDetectedFragments;
    }
    catch (NoRuleException | RulesException e ) 
    {
      msnAnalyzer_ = new MSnAnalyzer(lipidClassName_, lipidAdduct_, data_, analyzer_, null,true,false);
      return 0;
    }      
  }
  
  /**
   * Checks the equation of the head tab before add
   * @throws RulesException 
   * @throws CgException 
   * @throws SpectrummillParserException 
   * @throws IOException 
   * @throws NoRuleException 
   * @throws ChemicalFormulaException thrown if there is something wrong with the formula
   */
  public void checkToAddHeadIntensityRule() throws RulesException, IOException, SpectrummillParserException, CgException, NoRuleException, ChemicalFormulaException
  {
    //if (!this.checkGeneralEntries(true)) return;
    int before = 0;
    int after = 0;

    //check if the equation field is not empty
    if (equRuleEquationField_.getText()==null || equRuleEquationField_.getText().length()==0) return;
    
    headRuleTab_.remove( equRuleOK_ );            
    equRuleOK_ = new JLabel( new ImageIcon( redCrossPicture_ ));   
    gridBayLayout_.setConstraints(equRuleOK_, new GridBagConstraints (0, headIntensityRules_.size()+headFragments_.size()+5, 1, 1, 0, 0, 
    GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0 ));
    headRuleTab_.add( equRuleOK_ );
    headRuleTab_.repaint();
    headRuleTab_.validate();  
    
    //The container has to be set back
    paintNewSpectra(false); 
    
    before = headIntensityCounter();    
    IntensityRuleVO ruleVO;     
    ruleVO = FragRuleParser.extractIntensityVOFromEquation(equRuleEquationField_.getText(), 0, 2, FragmentRuleVO.getStringKeyHash(headFragmentRules_), FragmentRuleVO.getStringKeyHash(chainFragmentRules_), generalSettingsVO_.getAmountOfChains());
    ruleVO.setMandatory(getIntensityMandatory(equRuleMandatoryCombo_.getSelectedIndex()));
    headIntensityRules_.add(ruleVO);    
    paintNewSpectra(false);
    after = headIntensityCounter();
    if(after ==  before+1)
    {       
      headRuleTab_.remove( equRuleOK_ );            
      equRuleOK_ = new JLabel( new ImageIcon( resource1_ ));   
      gridBayLayout_.setConstraints(equRuleOK_, new GridBagConstraints (0, headIntensityRules_.size()+headFragments_.size()+4, 1, 1, 0, 0, 
      GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0 ));
      headRuleTab_.add( equRuleOK_ );
      headRuleTab_.repaint();
      headRuleTab_.validate();  
    }
    
    headIntensityRules_.remove(ruleVO); 
  }
  
  /**
   * Counts the number of fulfilled head intensity rules
   * @return the number of the fulfilled intensity rules
   * @throws RulesException
   * @throws IOException
   * @throws SpectrummillParserException
   * @throws CgException when there is something wrong with quantitation
   * @throws ChemicalFormulaException thrown if there is something wrong with the formula
   */
  public int headIntensityCounter() throws RulesException, IOException, SpectrummillParserException, CgException, ChemicalFormulaException
  {
    int numberOfFulfilledIntensityRules = 0;   
    try 
    {
      msnAnalyzer_ = updateMSnAnalyzerToCurrentSettings();
    }
    catch (NoRuleException | HydroxylationEncodingException | LipidCombinameEncodingException e) 
    {
      e.printStackTrace();
    }  
    Hashtable<String,IntensityRuleVO> fulfilledIntensityRules = msnAnalyzer_.getFulfilledHeadIntensityRules();
    numberOfFulfilledIntensityRules = fulfilledIntensityRules.size();
    return numberOfFulfilledIntensityRules;
  }
  
  
  /**
   * Checks a new chain fragment
   * @throws RulesException
   * @throws IOException
   * @throws SpectrummillParserException
   * @throws CgException
   * @throws NoRuleException
   * @throws HydroxylationEncodingException thrown if hydroxylation the encoding does not exist
   * @throws ChemicalFormulaException thrown if there is something wrong with the formula
   * @throws LipidCombinameEncodingException thrown when a lipid combi id (containing type and OH number) cannot be decoded
   */
  public void checkToAddChainFragment() throws RulesException, IOException, SpectrummillParserException, CgException, NoRuleException, HydroxylationEncodingException, ChemicalFormulaException, LipidCombinameEncodingException 
  {
    //if (!this.checkGeneralEntries(true)) return;
    int before = 0;
    int after = 0;

    //check if the equation field is not empty
    if (chainRuleNameField_.getText()==null || chainRuleNameField_.getText().length()==0) return;
    
    this.chainTab_.remove( this.chainRuleOK_ );            
    this.chainRuleOK_ = new JLabel( new ImageIcon( this.redCrossPicture_ ));   
    this.gridBayLayout_.setConstraints(this.chainRuleOK_, new GridBagConstraints (0, this.chainFragmentRules_.size()+1, 1, 1, 0, 0, 
    GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0 ));
    this.chainTab_.add( this.chainRuleOK_ );
    this.chainTab_.repaint();
    this.chainTab_.validate();  
    
    //The container has to be set back
    paintNewSpectra(false); 
    
    before = chainFragmentCounter();
    FragmentRuleVO ruleVO = null;
    if(chainFragmentRules_.containsKey(chainRuleNameField_.getText())) return;
    if(checkChainFragment())
	    {
	      ruleVO = new FragmentRuleVO(this.chainRuleNameField_.getText(),this.chainRuleFormulaField_.getText(),Integer.parseInt(this.chainRuleChargeField_.getText()), 
	          Integer.parseInt(this.chainRuleMsLevelField_.getText()), getFragMandatory(this.chainRuleMandatoryCombo_.getSelectedIndex()), this.headFragmentRules_, this.chainFragmentRules_, this.elementParser_);
	      this.chainFragmentRules_.put(this.chainRuleNameField_.getText(), ruleVO); 
	      this.sortedChainFragments_ = getSortedVector(this.chainFragmentRules_, TYPE_CHAIN); 
	    } else return;
	    paintNewSpectra(false); 
	    after = chainFragmentCounter();
	    if(after ==  before+1)
	    {    
	      this.chainTab_.remove( this.chainRuleOK_ ); 
	      this.chainRuleOK_ = new JLabel( new ImageIcon( this.resource1_ ));  
	      this.gridBayLayout_.setConstraints(this.chainRuleOK_, new GridBagConstraints (0, this.chainFragmentRules_.size(), 1, 1, 0, 0, 
	      GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0 ));
	      this.chainTab_.add( chainRuleOK_ );
	      this.chainTab_.repaint();
	      this.chainTab_.validate();  
	    }  
	    
	    if(checkChainFragment())
	    {
	    	this.chainFragmentRules_.remove(this.chainRuleNameField_.getText());
	      this.sortedChainFragments_ = getSortedVector(this.chainFragmentRules_, TYPE_CHAIN);
	    }
  }
  
  /**
   * Checks if there is any problem while building the chain fragmentVO
   * @return true for no problems, false for an error happened
   */
  public boolean checkChainFragment()
  {
    boolean check = true;
    try 
    {
      new FragmentRuleVO(chainRuleNameField_.getText(),chainRuleFormulaField_.getText(),Integer.parseInt(chainRuleChargeField_.getText()), 
      Integer.parseInt(chainRuleMsLevelField_.getText()), getFragMandatory(chainRuleMandatoryCombo_.getSelectedIndex()), headFragmentRules_, chainFragmentRules_, elementParser_);
    }
    catch (NumberFormatException e) 
    {
      check = false;
    }
    catch (RulesException e) 
    {
      check = false;
    }
    if (chainRuleNameField_.getText()==null || chainRuleNameField_.getText().length()==0 ||
        chainRuleFormulaField_.getText()==null || chainRuleFormulaField_.getText().length()==0 ||
        chainRuleChargeField_.getText()==null || chainRuleChargeField_.getText().length()==0 ||
        chainRuleMsLevelField_.getText()==null || chainRuleMsLevelField_.getText().length()==0)
      check = false;
    if (!FragmentVerifier.checkFragmentNameInput(chainRuleNameField_)) check = false;
    if (!FragmentVerifier.checkFragmentFormulaInput(chainRuleFormulaField_)) check = false;
      
    return check;
  }
  
  /**
   * Counts the maximum of detected chain fragments
   * @return the maximum of detected chain fragments
   * @throws RulesException
   * @throws IOException
   * @throws SpectrummillParserException
   * @throws CgException
   * @throws HydroxylationEncodingException thrown if hydroxylation the encoding does not exist
   * @throws ChemicalFormulaException thrown if there is something wrong with the formula
   * @throws LipidCombinameEncodingException  
   */
  public int chainFragmentCounter() throws RulesException, IOException, SpectrummillParserException, CgException, HydroxylationEncodingException, ChemicalFormulaException, LipidCombinameEncodingException 
  {
    int maxFragmentsDetected = 0;    
    try 
    {
      msnAnalyzer_ = this.updateMSnAnalyzerToCurrentSettings();
      Hashtable<String,Hashtable<String,CgProbe>> chainFragments = msnAnalyzer_.getChainFragments();  
      Enumeration<String> chainEnum = chainFragments.keys(); 
      while (chainEnum.hasMoreElements()) 
      {        
        String currentKey = (String)chainEnum.nextElement();      
        Hashtable<String,CgProbe> chainFragmentsDetected = (Hashtable<String,CgProbe>) chainFragments.get(currentKey);
        if(maxFragmentsDetected < chainFragmentsDetected.size())
        {
          maxFragmentsDetected = chainFragmentsDetected.size();
        }     
      }

    }
    catch (NoRuleException | RulesException e) 
    {
      msnAnalyzer_ = new MSnAnalyzer(lipidClassName_, lipidAdduct_, data_, analyzer_, null, true, false);
      return 0;
    } 
    
    return maxFragmentsDetected;
  }
  

  /**
   * Checks a new chain equation if it was found/used before it puts it to the others
   * @throws RulesException
   * @throws IOException
   * @throws SpectrummillParserException
   * @throws NoRuleException
   * @throws CgException
   * @throws ChemicalFormulaException thrown if there is something wrong with the formula
   */
  public void checkToAddChainIntensityRule() throws RulesException, IOException, SpectrummillParserException, NoRuleException, CgException, ChemicalFormulaException
  {    
    //if (!this.checkGeneralEntries(true)) return;
    int before = 0;
    int after = 0;

    chainTab_.remove( equChainRuleOK_ );            
    equChainRuleOK_ = new JLabel( new ImageIcon( redCrossPicture_ ));   
    gridBayLayout_.setConstraints(equChainRuleOK_, new GridBagConstraints (0, chainIntensityRules_.size()+chainFragmentRules_.size()+5, 1, 1, 0, 0, 
    GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0 ));
    chainTab_.add( equChainRuleOK_ );
    chainTab_.repaint();
    chainTab_.validate();  
    //The container has to be set back
    paintNewSpectra(false); 
    
    before = chainIntensityCounter();    
    IntensityRuleVO ruleVO;     
    ruleVO = FragRuleParser.extractIntensityVOFromEquation(equChainRuleEquationField_.getText(), 0, 2, FragmentRuleVO.getStringKeyHash(headFragmentRules_), FragmentRuleVO.getStringKeyHash(chainFragmentRules_), generalSettingsVO_.getAmountOfChains());
    ruleVO.setMandatory(getIntensityMandatory(equChainRuleMandatoryCombo_.getSelectedIndex()));
    chainIntensityRules_.add(ruleVO);    
    
    
    paintNewSpectra(false);
    after = chainIntensityCounter();
    if(after ==  before+1)
    {      
      chainTab_.remove( equChainRuleOK_ );            
      equChainRuleOK_ = new JLabel( new ImageIcon( resource1_ ));   
      gridBayLayout_.setConstraints(equChainRuleOK_, new GridBagConstraints (0, chainIntensityRules_.size()+chainFragmentRules_.size()+4, 1, 1, 0, 0, 
      GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0 ));
      chainTab_.add( equChainRuleOK_ );
      chainTab_.repaint();
      chainTab_.validate();  
    }
    chainIntensityRules_.remove(ruleVO); 
  }
  
  /**
   * Counts the number of fulfilled chain intensity rules
   * @return
   * @throws RulesException
   * @throws IOException
   * @throws SpectrummillParserException
   * @throws CgException
   * @throws ChemicalFormulaException thrown if there is something wrong with the formula
   */
  public int chainIntensityCounter() throws RulesException, IOException, SpectrummillParserException, CgException, ChemicalFormulaException
  {    
    try 
    {
      msnAnalyzer_ = this.updateMSnAnalyzerToCurrentSettings();
    }
    catch (NoRuleException | HydroxylationEncodingException | LipidCombinameEncodingException e) 
    {
      e.printStackTrace();
    }   
    Hashtable<String,Hashtable<String,IntensityChainVO>> fulfilledChainIntensityRules = msnAnalyzer_.getFulfilledChainIntensityRules();
    
    int counter = 0;
    Enumeration<String> chainEnum = fulfilledChainIntensityRules.keys();         
    while (chainEnum.hasMoreElements()) 
    {        
      String currentKey = (String)chainEnum.nextElement();      
      Hashtable<String,IntensityChainVO> currentfulfilledChainIntensityRules = fulfilledChainIntensityRules.get(currentKey);
      if(currentfulfilledChainIntensityRules.size() > counter)
        counter = currentfulfilledChainIntensityRules.size();      
    }    
    return counter;  
  }
  

  /**
   * Checks a new equation if it was found/used before it puts it to the others
   * @throws SpectrummillParserException 
   * @throws IOException 
   * @throws RulesException 
   * @throws NoRuleException 
   * @throws CgException 
   */
  public void checkToAddPositionRule() throws RulesException, IOException, SpectrummillParserException, NoRuleException, CgException
  {
    //if (!this.checkGeneralEntries(true)) return;
    int before = 0;
    int after = 0;

    //check if the equation field is not empty
    if (positionEquationField_.getText()==null || positionEquationField_.getText().length()==0) return;
    
    positionTab_.remove( positionRuleOK_ );            
    positionRuleOK_ = new JLabel( new ImageIcon( redCrossPicture_ ));     
    gridBayLayout_.setConstraints(positionRuleOK_, new GridBagConstraints (0, 1+positionIntensityRules_.size(), 1, 1, 0, 0, 
    GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0 ));
    positionTab_.add( positionRuleOK_ );
    positionTab_.repaint();
    positionTab_.validate();  
    //The container has to be set back
    paintNewSpectra(false); 
    
    before = unfulfilledPosRulesCounter();
    try
    {
      IntensityRuleVO ruleVO = FragRuleParser.extractIntensityVOFromEquation(positionEquationField_.getText(), 0, FragRuleParser.POSITION_SECTION, FragmentRuleVO.getStringKeyHash(headFragmentRules_), FragmentRuleVO.getStringKeyHash(chainFragmentRules_),generalSettingsVO_.getAmountOfChains());         
      ruleVO.setMandatory(getIntensityMandatory(positionRuleMandatoryCombo_.getSelectedIndex()));
      positionIntensityRules_.add(ruleVO);     
      saveRules(CACHE_DIR);
      after = unfulfilledPosRulesCounter();    
      if(!(positionEquationField_.getText().equals("")))
      {
        if(after ==  before)
        {   
                  
          positionTab_.remove( positionRuleOK_ );            
          positionRuleOK_ = new JLabel( new ImageIcon( resource1_ ));     
          gridBayLayout_.setConstraints(positionRuleOK_, new GridBagConstraints (0, positionIntensityRules_.size(), 1, 1, 0, 0, 
              GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0 ));
          positionTab_.add( positionRuleOK_ );
          positionTab_.repaint();
          positionTab_.validate();   
        }       
        positionIntensityRules_.remove(ruleVO);       
      }
      else
      {
        positionIntensityRules_.remove(ruleVO);
      }
    }
    catch (RulesException e2) 
    { 
    } 
  }
  
  /**
   * Counts the unfulfilled position rules
   * @return the number of unfulfilled position rules
   */
  public int unfulfilledPosRulesCounter() 
  { 
    int result = 0;
    try 
    {
      MSnDebugVO debugInfo;
      debugInfo = getCurrentDebugInfo();
      Hashtable<String,Hashtable<String,IntensityRuleVO>> unfulfilledPosRules = debugInfo.getUnfulfilledPositionRules();
      result = unfulfilledPosRules.size();
    }
    catch (RulesException | IOException | SpectrummillParserException | CgException | HydroxylationEncodingException
        | ChemicalFormulaException | LipidCombinameEncodingException e) 
    {
      e.printStackTrace();
    }
    return result;    
  }
  
  
  /**
   * Sets the input field to null to create a new class
   */
  private void setEverythingNull()
  {    
    generalSettingsVO_ = new GeneralSettingsVO();
  
    headFragmentRules_= new Hashtable<String, FragmentRuleVO>();                 
    headIntensityRules_ = new Vector<IntensityRuleVO>();
    chainFragmentRules_= new Hashtable<String, FragmentRuleVO>();        
    chainIntensityRules_ = new Vector<IntensityRuleVO>();
   positionIntensityRules_ = new Vector<IntensityRuleVO>();
  
    if (topSplitPane_!=null){
      topSplitPane_.remove( middleTabSection_ );
      topSplitPane_.validate();
      ruleTabsSection();
      middleTabSection_.setSelectedIndex(0);
      topSplitPane_.setBottomComponent(middleTabSection_);
    }
  }  
  
  /**
   * Returns a sorted string of results
   * @param result
   * @param showRetentionTime
   * @return the sorted result string
   */
  public String getResultWriteString(LipidParameterSet result, boolean showRetentionTime)
  {
//    try {
//      msnAnalyzer_ = new MSnAnalyzer(lipidClassName_, lipidAdduct_, data_,
//          analyzer_, null, true, false);
//    }
//    catch (RulesException | IOException | SpectrummillParserException
//        | CgException | HydroxylationEncodingException | ChemicalFormulaException | LipidCombinameEncodingException e) {
//      e.printStackTrace();
//    }
//
//    LipidParameterSet result = msnAnalyzer_.getResult();

    String resultString = "";
    LipidomicsMSnSet msnSet = (LipidomicsMSnSet) result;
    Vector<String> detected = msnSet.getMSnIdentificationNamesWithSNPositions();
    for (String name : detected) {
        resultString = resultString + name + "_" + result.getRt()
            + "<br>";
    }
    return resultString;
  }

  /**
   * checks if the fragment is OK before updating the hash tables
   * @param ruleName
   * @param formula
   * @param charge
   * @param msLevel
   * @param mandatory
   * @param type 1...Formula, 2...Charge, 3... MSLevel, 4... Mandatory
   * @param fragmentType 1...HeadFragments, 2...ChainFragments
   * @throws RulesException
   */
  public FragmentRuleVO checkFragment(String ruleName, String formula, int charge, int msLevel, short mandatory, int type, int fragmentType) throws RulesException
  {
    FragmentRuleVO ruleVO = null;
    Hashtable<String,FragmentRuleVO> source = null;
    if(fragmentType == 1)
    {
      source = headFragmentRules_;
    }
    if(fragmentType == 2)
    {
      source = chainFragmentRules_;
    }
    
    Enumeration<String> headFragmentEnum = source.keys(); 
    while (headFragmentEnum.hasMoreElements()) 
    {
      
      String currentKey = (String)headFragmentEnum.nextElement();      
      FragmentRuleVO rule = (FragmentRuleVO) source.get(currentKey);  
     
        if(currentKey.equals(ruleName))
        {        
          if(type == 1)
          {
            ruleVO = new FragmentRuleVO(currentKey,formula,rule.getCharge(), rule.getMsLevel(), mandatory, headFragmentRules_, chainFragmentRules_, elementParser_);          
          }
          if(type == 2)
          {
            ruleVO = new FragmentRuleVO(currentKey,formula,charge, rule.getMsLevel(), mandatory, headFragmentRules_, chainFragmentRules_, elementParser_);          
          }
          if(type == 3)
          {
            ruleVO = new FragmentRuleVO(currentKey, formula, charge , msLevel, mandatory, headFragmentRules_, chainFragmentRules_, elementParser_);          
          }         
          if (type==1 || type==2 || type==3)
            return ruleVO;
	        }
    } 
    return null;
  }
  
  
  /**
   * updates the hash tables and the sorted vectors with the changed FragementRuleVO
   * @param ruleVO the ruleVO to update
   * @param fragmentType the type (TYPE_HEAD or TYPE_CHAIN)
   */
  public void updateRuleVO(FragmentRuleVO ruleVO, int fragmentType){
    Hashtable<String,FragmentRuleVO> source = null;
    if(fragmentType == 1)
    {
      source = headFragmentRules_;
    }
    if(fragmentType == 2)
    {
      source = chainFragmentRules_;
    }
    source.put(ruleVO.getName(),ruleVO);
    headFragments_ = getSortedVector(headFragmentRules_, TYPE_HEAD);
    sortedChainFragments_ = getSortedVector(chainFragmentRules_, TYPE_CHAIN);
  }
  
  /**
   * checks if there are any fulfilled Spectra with the currently entered rules
   * @param specNumber the number of the selected spectrum
   * @return the identified MSn fragments
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   * @throws CgException errors from the quantitation process
   * @throws HydroxylationEncodingException HydroxylationEncodingException
   * @throws ChemicalFormulaException thrown if there is something wrong with the formula
   * @throws LipidCombinameEncodingException thrown when a lipid combi id (containing type and OH number) cannot be decoded
   */
  public LipidParameterSet testForMSnDetection(int specNumber) throws RulesException, NoRuleException, IOException, SpectrummillParserException, CgException, HydroxylationEncodingException, ChemicalFormulaException, LipidCombinameEncodingException{
    try{
      MSnAnalyzer analyzer = updateMSnAnalyzerToCurrentSettings(specNumber);
      data_ = analyzer.getResult();
    } catch (RulesException rx){
      if (rx.getMessage().endsWith(FragRuleParser.NO_HEAD_AND_CHAINS_SECTION) || !this.showFragmentTabs_){
        data_ = new LipidParameterSet(data_);
      } else throw rx;
    }
    return data_;
  }
  
  /**
   * creates an MSnAnalyzer object with the currently entered rules
   * @return the MSnAnalyzer with the currently entered rules
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   * @throws CgException errors from the quantitation process
   * @throws HydroxylationEncodingException thrown if hydroxylation the encoding does not exist
   * @throws ChemicalFormulaException thrown if there is something wrong with the formula
   * @throws LipidCombinameEncodingException thrown when a lipid combi id (containing type and OH number) cannot be decoded
   */
  private MSnAnalyzer updateMSnAnalyzerToCurrentSettings() throws RulesException, NoRuleException, IOException, SpectrummillParserException, CgException, HydroxylationEncodingException, ChemicalFormulaException, LipidCombinameEncodingException{
    return this.updateMSnAnalyzerToCurrentSettings(spectrumUpdater_.getMs2LevelSpectrumSelected());
  }
  
  /**
   * creates an MSnAnalyzer object with the currently entered rules for a selected spectrum
   * @param specNumber the number of the selected spectrum
   * @return the MSnAnalyzer with the currently entered rules
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws NoRuleException thrown if the rules are not there
   * @throws IOException exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   * @throws CgException errors from the quantitation process
   * @throws HydroxylationEncodingException thrown if hydroxylation the encoding does not exist
   * @throws ChemicalFormulaException thrown if there is something wrong with the formula
   * @throws LipidCombinameEncodingException thrown when a lipid combi id (containing type and OH number) cannot be decoded
   */
  private MSnAnalyzer updateMSnAnalyzerToCurrentSettings(int specNumber) throws RulesException, NoRuleException, IOException, SpectrummillParserException, CgException, HydroxylationEncodingException, ChemicalFormulaException, LipidCombinameEncodingException{
   FragmentCalculator fragCalc_ = new FragmentCalculator(CACHE_DIR,lipidClassName_,lipidAdduct_,data_.getNameStringWithoutRt(),data_.getChemicalFormula(),
        data_.getChemicalFormulaWODeducts(),data_.Mz[0],data_.getCharge(),data_.getOhNumber(),data_.getOxState());
   	float tol = LipidomicsConstants.getMs2PrecursorTolerance(data_.Mz[0]);
    analyzer_.prepareMSnSpectraCache(data_.Mz[0]-tol, data_.Mz[0]+tol, LipidomicsConstants.getMs2MinIntsForNoiseRemoval());
    Vector<Range> ranges = analyzer_.findSingleSpectraRanges(fragCalc_.getSpectrumLevelRange());     
    Vector<CgProbe> probes = new Vector<CgProbe>();
    if (specNumber<0 || specNumber>=ranges.size()){
      for (Range range : ranges)
      {
        CgProbe probe = createDummyProbe(data_.Mz[0]);
        probe.Peak = (range.getStart()+range.getStop())/2f;
        probe.LowerValley = range.getStart();
        probe.UpperValley = range.getStop();
        probes.add(probe);
      }
    } else {
      Range range = ranges.get(specNumber);
      CgProbe probe = createDummyProbe(data_.Mz[0]);
      probe.Peak = (range.getStart()+range.getStop())/2f;
      probe.LowerValley = range.getStart();      
      probe.UpperValley = range.getStop();      
      probes.add(probe);
    }
    data_.setProbes(probes);
    MSnAnalyzer analyzer = new MSnAnalyzer(CACHE_DIR, lipidClassName_, lipidAdduct_, data_, analyzer_, null, false, true, true);
    return analyzer;
  }
  
  /**
   * creates an empty CgProbe object for defining retention time ranges
   * @param mz the m/z value of the CgProbe object
   * @return the empty CgProbe object for defining retention time ranges
   */
  private CgProbe createDummyProbe(float mz){
    CgProbe probe = new CgProbe(0,1);
    probe.AreaStatus = CgAreaStatus.OK;
    probe.Area = 100f;
    probe.AreaError = 0f;
    probe.Background = 0f;
    probe.Mz = mz;
    probe.LowerMzBand = LipidomicsConstants.getCoarseChromMzTolerance(mz);
    probe.UpperMzBand = LipidomicsConstants.getCoarseChromMzTolerance(mz);
    probe.isotopeNumber = 0;
    return probe;
  }
 
  /**
   * Paints a new Spectrum
   * @throws RulesException 
   */
  public void paintNewSpectra(boolean newPrec)
  {      
    try 
    {
      saveRules(CACHE_DIR);
      spectrumUpdater_.updateSpectra(data_, newPrec);
//      lipidPrecursorValueBefore_ = lipidPrecursorValue_;
//      lipidNameBefore_ = lipidName_;
    }
    catch (CgException | IOException | RulesException | SpectrummillParserException
        | NoRuleException e) 
    {     
      e.printStackTrace();
    }
  }  
  
  

  /**
   * Adds Equation for the verifier
   * @param ruleVO the IntensityRuleVO
   * @param type 1... head, 2... chain, 3... position
   * @param position the position in the vector
   */
  public void addEquationToIntensityRules(IntensityRuleVO ruleVO, int type, int position)
  {
    if(type == 1)
      headIntensityRules_.set(position, ruleVO);  
    
    if(type == 2)
      chainIntensityRules_.set(position, ruleVO);
    
    if(type == 3)
      positionIntensityRules_.set(position, ruleVO);
    
  }
  
  /**
   * returns the head rules for the verifier
   */
  public  Hashtable<String,FragmentRuleVO> getHeadFragmentRules()
  {
    return headFragmentRules_;
  }
  
  /**
   * returns the chain fragments for the verifier
   * @return
   */
  public  Hashtable<String,FragmentRuleVO> getChainFragmentRules()
  {
    return chainFragmentRules_;
  }
  
  /**
   * Returns the head fragments sorted vector
   * @return
   */
  public Vector<FragmentRuleVO> getHeadFragmentVector()
  {
	  return headFragments_;
  }
  
  /**
   * Returns the sorted chain fragments in a vector
   * @return
   */
  public Vector<FragmentRuleVO> getChainFragmentVector()
  {
	  return sortedChainFragments_;
  }
    
  /**
   * Return the amount of chains for the verifier
   * @return
   */
  public String getAmountOfChains()
  {
    return generalSettingsVO_.getAmountOfChains().toString();
  }
  
  public boolean intToBoolean(int input)
  {
  	if(input == 0)
  	{
  		return true;
  	}  		
  	if(input == 1)
  	{
  		return false;
  	}  		
		return false;
  }
  
  /**
   * Returns the head formula fields
   * @return
   */
  public JTextField[] getHeadRuleFormulas()
  {
	  return headRuleFormulas_;
  }   
 
  /**
   * Returns the head charge fields
   * @return
   */
  public JTextField[] getHeadRuleCharges()
  {
	  return headRuleCharges_;
  }  
  
  /**
   * Returns the head mslevel fields
   * @return
   */
  public JTextField[] getHeadRuleMSLevels()
  {
	  return headRuleMSLevels_;
  }  
  
  /**
   * Returns the chain rule formula fields
   * @return
   */
  public JTextField[] getChainRuleFormulas()
  {
	  return chainRuleFormulas_;
  }   
 
  /**
   * Returns the chain charge fields
   * @return
   */
  public JTextField[] getChainRuleCharges()
  {
	  return chainRuleCharges_;
  }  
  
  /**
   * Returns the chain rule mslevel fields
   * @return
   */
  public JTextField[] getChainRuleMSLevels()
  {
	  return chainRuleMSLevels_;
  }  
  
  /**
   * If the focus is lost it checks if the entry is correct
   * @param position
   * @param type
   * @param headOrChain
   */
  public void checkFragmentWithErrors(int position, int type, int headOrChain)
  {
//    if (!this.checkGeneralEntries()) return;
    if(headOrChain == 1)
  	{
		  if(type == 0)
		  {
			  headFragmentFormulaVerifier_[position].checkFragment();	  
		  }
		  
		  if(type == 1)
		  {
			  headFragmentChargeVerifier_[position].checkFragment();	 
		  }
		  
		  if(type == 2)
		  {
			  headFragmentMsLevelVerifier_[position].checkFragment();	 
		  }
  	}
  	if(headOrChain == 2)
  	{
		  if(type == 0)
		  {
			  chainFragmentFormulaVerifier_[position].checkFragment();	  
		  }
		  
		  if(type == 1)
		  {
			  chainFragmentChargeVerifier_[position].checkFragment();	 
		  }
		  
		  if(type == 2)
		  {
			  chainFragmentMsLevelVerifier_[position].checkFragment();	 
		  }
  	}
  }

  /**
   * While typing it checks if it is correct and saves it in the caches
   * @param position
   * @param type
   * @param headOrChain
   */
  public void refreshTextfieldToVerifie(int position, int type, int headOrChain)
  {
    //if (!this.checkGeneralEntries(true)) return;
  if(headOrChain == 1)
  	{
		  if(type == 0)
		  {
			  headFragmentFormulaVerifier_[position].refresh(position);	  
		  }
		  
		  if(type == 1)
		  {
			  headFragmentChargeVerifier_[position].refresh(position);	 
		  }
		  
		  if(type == 2)
		  {
			  headFragmentMsLevelVerifier_[position].refresh(position);	 
		  }
  	}
  	if(headOrChain == 2)
  	{
		  if(type == 0)
		  {
			  chainFragmentFormulaVerifier_[position].refresh(position);	  
		  }
		  
		  if(type == 1)
		  {
			  chainFragmentChargeVerifier_[position].refresh(position);	 
		  }
		  
		  if(type == 2)
		  {
			  chainFragmentMsLevelVerifier_[position].refresh(position);	 
		  }
  	}
  }
  
  /**
   * Refreshes the variables of the verifiers for the first time
   * @param position
   * @param type
   * @param headOrChain
   */
  public void refreshTextfieldToVerifieFirstTime(int position, int type, int headOrChain)
  {	  
  	if(headOrChain == 1)
  	{
  		 if(type == 0)
  		  {
  			  headFragmentFormulaVerifier_[position].refreshFirstTime(position);	  
  		  }
  		  
  		  if(type == 1)
  		  {
  			  headFragmentChargeVerifier_[position].refreshFirstTime(position);	 
  		  }
  		  
  		  if(type == 2)
  		  {
  			  headFragmentMsLevelVerifier_[position].refreshFirstTime(position);	 
  		  }
  	}
  	if(headOrChain == 2)
  	{
  		 if(type == 0)
  		  {
  			  chainFragmentFormulaVerifier_[position].refreshFirstTime(position);	  
  		  }
  		  
  		  if(type == 1)
  		  {
  			  chainFragmentChargeVerifier_[position].refreshFirstTime(position);	 
  		  }
  		  
  		  if(type == 2)
  		  {
  			  chainFragmentMsLevelVerifier_[position].refreshFirstTime(position);	 
  		  }
  	}
  } 
  
  /**
   * Refreshes the verifiers of the equations for the first time
   * @param position
   * @param headOrChainOrPostion
   */
  public void refreshEquationTextfieldToVerifieFirstTime(int position,  int headOrChainOrPostion)
  {	  
  	if(headOrChainOrPostion == 1)
  	{
  		headEquationVerifier_[position].refreshFirstTime();
  	}
  	if(headOrChainOrPostion == 2)
  	{
  		chainFragmentFormulaVerifier_[position].refreshFirstTime(position);
  	}
  	if(headOrChainOrPostion == 3)
  	{
  		chainFragmentFormulaVerifier_[position].refreshFirstTime(position);	
  	}
  }
  
  /**
   * Returns the head fragment mandatory values
   * @return
   */
  public boolean[] getHeadFragmentMandatories()
  {
  	return this.headFragmentMandatories_;
  }
  
  /**
   * Returns the chain fragment mandatory values
   * @return
   */
  public boolean[]   getChainFragmentMandatories()
  {
  	return this.chainFragmentMandatories_;
  }
  
  /**
   * Returns the head equations mandatory values
   * @return
   */
  public boolean[] getHeadEquationtMandatories()
  {
  	return this.headEquationMandatories_;
  }
  

  /**
   * Checks the equations while they are changed and saves it in the caches
   * @param position
   * @param headOrChainOrPosition
   */
	  public void refreshEquationFieldToVerifie(int position, int headOrChainOrPosition)
	  {
    //if (!this.checkGeneralEntries(true)) return;
  if(headOrChainOrPosition == 1)
		{
			headEquationVerifier_[position].refresh();
		}
		if(headOrChainOrPosition == 2)
		{
			chainEquationVerifier_[position].refresh();	
		}
		if(headOrChainOrPosition == 3)
		{
			positionEquationVerifier_[position].refresh();	
		}
	  } 
	
	/**
	 * Returns the head equation fields
	 * @return
	 */
	public JTextField[] getHeadEquations()
	{
		return this.headRuleEquationFields_;
	}
	
	/**
	 * Returns the chain equation fields
	 * @return
	 */
	public JTextField[] getChainEquations()
	{
		return this.chainRuleEquationFields_;
	}
	
	/**
	 * Returns the chain equations mandatory fields
	 * @return
	 */
	public boolean[] getChainEquationtMandatories()
	{
	 	return this.chainEquationMandatories_;
	}
	 
	/**
	 * Returns the position equations fields
	 * @return
	 */
	public JTextField[] getPositionEquations()
	{
		return this.positionRuleEquationFields_;
	}
		
	/**
	 * Returns the position equations mandatory fields
	 * @return
	 */
	public boolean[] getPositionEquationtMandatories()
	{
		return this.positionEquationMandatories_;
	}

  
  /**
   * Constructor of the rule definition class
   * @param lipidClassName the name of the lipid class
   * @param lipidName_ the name of the lipid
   * @param lipidPrecursorValue_ the current value of the precursor
   * @throws ChemicalFormulaException 
   * @throws thrown if the rules are not there
   */
  public RuleDefinitionInterface(String lipidClassName, LipidParameterSet data, LipidomicsAnalyzer analyzer, int highestMSLevel, SpectrumUpdateListener object) throws NoRuleException, ChemicalFormulaException 
  {    
    super(JSplitPane.VERTICAL_SPLIT);    
    showFragmentTabs_ = true;
    setDividerSize(0);    
    this.data_ = data;    
    this.lipidClassName_ = lipidClassName;    
    this.lipidName_ = data.getNameString();
    this.lipidPrecursorValue_ = data.Mz[0];
    this.lipidFormula_ = StaticUtils.getFormulaInHillNotation(StaticUtils.categorizeFormula(data.getAnalyteFormula()),true);
    this.lipidAdduct_ = data.getModificationName();    
    this.ruleClassIdentifier_ = StaticUtils.getRuleName(lipidClassName_, lipidAdduct_);
    this.analyzer_ = analyzer;
    this.spectrumUpdater_ = object;
//    this.lipidPrecursorValueBefore_ = lipidPrecursorValue_;
//    this.lipidNameBefore_ = lipidName_;    
    this.highestMSLevel_ = highestMSLevel;
    try
    {      
      RulesContainer.clearCache();
      getRulesFromCacheFirst();
    }    
    catch (SpectrummillParserException e) 
    {
      JOptionPane.showMessageDialog(topSplitPane_, e.getMessage()); 
    }
    catch (RulesException e) 
    {
      JOptionPane.showMessageDialog(topSplitPane_, e.getMessage()); 
    }
    catch (IOException e) 
    {
      JOptionPane.showMessageDialog(topSplitPane_, e.getMessage());               
    }   
      this.createPanel();
      try {
        saveRules(CACHE_DIR);
      }
      catch (IOException | RulesException | SpectrummillParserException
          | NoRuleException e) {
        // TODO Auto-generated catch block
        e.printStackTrace();
      }
      this.setVisible(true);
      this.invalidate();
    }
  
  /**
   * clears the files in the rules cache directory
   */
  public static void clearCacheDir(){
    File dir = new File(CACHE_DIR);
    if (!dir.exists()) dir.mkdir();
    File[] files = dir.listFiles();
    for (File file : files){
      file.delete();
    }
  }

  
  /**
   * changes the visibility level of the rule definition tabs
   * @param visible true if they shall be made visible
   */
  private void setFragmentTabsVisible(boolean visible){
    this.headRuleTab_.setVisible(visible);
    this.chainTab_.setVisible(visible);
    this.positionTab_.setVisible(visible);
  }
  
  public boolean performStorageOfInitialGeneralSettings(){
    try {
      saveRules(CACHE_DIR);
      this.showFragmentTabs_ = true;
      return true;
    }
    catch (IOException | RulesException | SpectrummillParserException
        | NoRuleException e1) {
      e1.printStackTrace();
      return false;
    }
  }
  
  public void refreshGeneralSettingsDependantDisplay(){
    setFragmentTabsVisible(true);
    refreshMiddle(1);
  }
  public String getAnalyteName(){
    return data_.getNameStringWithoutRt();
  }
  
  public String getLipidClassName(){
    return ruleClassIdentifier_;
  }

  public boolean areGeneralSettingsDependingFieldsVisible()
  {
    return this.showFragmentTabs_;
  }
  
  public void updateGeneralSettings(){
    generalSettingsVO_ = this.generalSettings_.getValues();
    paintNewSpectra(false); 
  }
}






