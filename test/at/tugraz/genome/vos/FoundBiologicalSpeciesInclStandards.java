package at.tugraz.genome.vos;

import java.util.LinkedHashMap;

public class FoundBiologicalSpeciesInclStandards extends FoundBiologicalSpecies
{

  public static LinkedHashMap<String,LinkedHashMap<String,ReferenceInfoVO>> getPISpeciesOrbitrap(){
    LinkedHashMap<String,LinkedHashMap<String,ReferenceInfoVO>> species = FoundBiologicalSpecies.getPISpeciesOrbitrap();
    LinkedHashMap<String,ReferenceInfoVO> molSpecies = new LinkedHashMap<String,ReferenceInfoVO>();
    molSpecies.put("12:0/13:0", new ReferenceInfoVO("12:0/13:0",new double[]{13.0d},true,true));
    species.put("25:0", molSpecies);
    molSpecies = new LinkedHashMap<String,ReferenceInfoVO>();
    molSpecies.put("17:0/14:1", new ReferenceInfoVO("17:0/14:1",new double[]{20.1d},true,true));
    species.put("31:1", molSpecies);
    molSpecies = new LinkedHashMap<String,ReferenceInfoVO>();
    molSpecies.put("21:0/22:6", new ReferenceInfoVO("21:0/22:6",new double[]{27.9d},true,true));
    species.put("43:6", molSpecies);
    return species;
  }
  
  public static LinkedHashMap<String,LinkedHashMap<String,ReferenceInfoVO>> getPPCSpeciesOrbitrap(){
    LinkedHashMap<String,LinkedHashMap<String,ReferenceInfoVO>> species = new LinkedHashMap<String,LinkedHashMap<String,ReferenceInfoVO>>();
    LinkedHashMap<String,ReferenceInfoVO> molSpecies = new LinkedHashMap<String,ReferenceInfoVO>();
    molSpecies.put("P-18:0/18:1", new ReferenceInfoVO("P-18:0/18:1",new double[]{29.8d},true,true));
    species.put("36:1", molSpecies);
    molSpecies = new LinkedHashMap<String,ReferenceInfoVO>();
    molSpecies.put("P-18:0/20:4", new ReferenceInfoVO("P-18:0/20:4",new double[]{28.0d},true,true));
    species.put("38:4", molSpecies);
    molSpecies = new LinkedHashMap<String,ReferenceInfoVO>();
    molSpecies.put("P-18:0/22:6", new ReferenceInfoVO("P-18:0/22:6",new double[]{27.0d},true,true));
    species.put("40:6", molSpecies);
    return species;
  }
  
  public static LinkedHashMap<String,LinkedHashMap<String,ReferenceInfoVO>> getPPESpeciesOrbitrap(){
    return FoundBiologicalSpecies.getPPESpeciesOrbitrap();
  }
  
  public static LinkedHashMap<String,LinkedHashMap<String,ReferenceInfoVO>> getLPCSpeciesOrbitrap(){
    LinkedHashMap<String,LinkedHashMap<String,ReferenceInfoVO>> species = FoundBiologicalSpecies.getLPCSpeciesOrbitrap();
    LinkedHashMap<String,ReferenceInfoVO> molSpecies = new LinkedHashMap<String,ReferenceInfoVO>();
    molSpecies.put("13:0", new ReferenceInfoVO("13:0",new double[]{1.7d},true,false));
    species.put("13:0", molSpecies);
    molSpecies = new LinkedHashMap<String,ReferenceInfoVO>();
    molSpecies.put("14:0", new ReferenceInfoVO("14:0",new double[]{2.4d},true,false));
    species.put("14:0", molSpecies);
    return species;
  }
  
  public static LinkedHashMap<String,LinkedHashMap<String,ReferenceInfoVO>> getLPESpeciesOrbitrap(){
    LinkedHashMap<String,LinkedHashMap<String,ReferenceInfoVO>> species = FoundBiologicalSpecies.getLPESpeciesOrbitrap();
    LinkedHashMap<String,ReferenceInfoVO> molSpecies = new LinkedHashMap<String,ReferenceInfoVO>();
    molSpecies.put("14:0", new ReferenceInfoVO("14:0",new double[]{2.5d},true,false));
    species.put("14:0", molSpecies);    
    return species;
  }
  
  public static LinkedHashMap<String,LinkedHashMap<String,ReferenceInfoVO>> getPSSpeciesOrbitrap(){
    LinkedHashMap<String,LinkedHashMap<String,ReferenceInfoVO>> species = FoundBiologicalSpecies.getPSSpeciesOrbitrap();
    LinkedHashMap<String,ReferenceInfoVO> molSpecies = new LinkedHashMap<String,ReferenceInfoVO>();
    molSpecies.put("16:0/16:0", new ReferenceInfoVO("16:0/16:0",new double[]{25.3d},true,false));
    species.put("32:0", molSpecies);
    molSpecies = new LinkedHashMap<String,ReferenceInfoVO>();
    molSpecies.put("17:0/17:0", new ReferenceInfoVO("17:0/17:0",new double[]{28.0d},true,false));
    species.put("34:0", molSpecies);
    molSpecies = new LinkedHashMap<String,ReferenceInfoVO>();
    molSpecies.put("16:0/18:1", new ReferenceInfoVO("16:0/18:1",new double[]{25.9d},true,true));
    species.put("34:1", molSpecies);
    molSpecies = new LinkedHashMap<String,ReferenceInfoVO>();
    molSpecies.put("16:0/18:2", new ReferenceInfoVO("16:0/18:2",new double[]{23.5d},true,true));
    species.put("34:2", molSpecies);
    molSpecies = new LinkedHashMap<String,ReferenceInfoVO>();
    molSpecies.put("18:0/18:0", new ReferenceInfoVO("18:0/18:0",new double[]{31.6d},true,false));
    species.put("36:0", molSpecies);
    molSpecies = new LinkedHashMap<String,ReferenceInfoVO>();
    molSpecies.put("18:0/18:1", new ReferenceInfoVO("18:0/18:1",new double[]{30.8d},true,true));
    species.put("36:1", molSpecies);    
    return species;
  }
  
  public static LinkedHashMap<String,LinkedHashMap<String,ReferenceInfoVO>> getLPSSpeciesOrbitrap(){
    LinkedHashMap<String,LinkedHashMap<String,ReferenceInfoVO>> species = new LinkedHashMap<String,LinkedHashMap<String,ReferenceInfoVO>>();
    LinkedHashMap<String,ReferenceInfoVO> molSpecies = new LinkedHashMap<String,ReferenceInfoVO>();
    molSpecies.put("16:0", new ReferenceInfoVO("16:0",new double[]{3.8d},true,false));
    species.put("16:0", molSpecies);
    molSpecies = new LinkedHashMap<String,ReferenceInfoVO>();
    molSpecies.put("17:1", new ReferenceInfoVO("17:1",new double[]{3.4d},true,false));
    species.put("17:1", molSpecies);
    molSpecies = new LinkedHashMap<String,ReferenceInfoVO>();
    molSpecies.put("18:0", new ReferenceInfoVO("18:0",new double[]{6.6d},true,false));
    species.put("18:0", molSpecies);
    molSpecies = new LinkedHashMap<String,ReferenceInfoVO>();
    molSpecies.put("18:1", new ReferenceInfoVO("18:1",new double[]{4.5d},true,false));
    species.put("18:1", molSpecies);
    return species;
  }
  
  public static LinkedHashMap<String,LinkedHashMap<String,ReferenceInfoVO>> getPCSpeciesOrbitrap(){
    LinkedHashMap<String,LinkedHashMap<String,ReferenceInfoVO>> species = FoundBiologicalSpecies.getPCSpeciesOrbitrap();
    species.get("32:0").put("14:0/18:0", new ReferenceInfoVO("14:0/18:0",new double[]{25.7d},true,true));
    LinkedHashMap<String,ReferenceInfoVO> molSpecies = new LinkedHashMap<String,ReferenceInfoVO>();
    molSpecies.put("20:0/20:0", new ReferenceInfoVO("20:0/20:0",new double[]{34.2d},true,false));
    species.put("40:0", molSpecies);
    molSpecies = new LinkedHashMap<String,ReferenceInfoVO>();
    molSpecies.put("24:1/24:1", new ReferenceInfoVO("24:1/24:1",new double[]{37.4d},true,false));
    species.put("48:2", molSpecies);
    return species;
  }
  
  public static LinkedHashMap<String,LinkedHashMap<String,ReferenceInfoVO>> getPESpeciesOrbitrap(){
    LinkedHashMap<String,LinkedHashMap<String,ReferenceInfoVO>> species = FoundBiologicalSpecies.getPESpeciesOrbitrap();
    LinkedHashMap<String,ReferenceInfoVO> molSpecies = new LinkedHashMap<String,ReferenceInfoVO>();
    molSpecies.put("16:0/16:0", new ReferenceInfoVO("16:0/16:0",26.1d,true,false));
    species.put("32:0", molSpecies);
    species.get("34:0").put("17:0/17:0", new ReferenceInfoVO("17:0/17:0",new double[]{28.6d},true,false));   
    molSpecies = new LinkedHashMap<String,ReferenceInfoVO>();
    molSpecies.put("18:0/18:0", new ReferenceInfoVO("18:0/18:0",30.8d,true,false));
    species.put("36:0", molSpecies);    
    return species;
  }
  

  public static LinkedHashMap<String,LinkedHashMap<String,ReferenceInfoVO>> getPGSpeciesOrbitrap(){
    LinkedHashMap<String,LinkedHashMap<String,ReferenceInfoVO>> species = FoundBiologicalSpecies.getPGSpeciesOrbitrap();
    LinkedHashMap<String,ReferenceInfoVO> molSpecies = new LinkedHashMap<String,ReferenceInfoVO>();
    molSpecies.put("16:0/16:0", new ReferenceInfoVO("16:0/16:0",24.2d,true,false));
    species.put("32:0", molSpecies);
    molSpecies = new LinkedHashMap<String,ReferenceInfoVO>();
    molSpecies.put("17:0/17:0", new ReferenceInfoVO("17:0/17:0",24.6d,true,false));
    species.put("34:0", molSpecies);
    molSpecies = new LinkedHashMap<String,ReferenceInfoVO>();
    molSpecies.put("18:0/18:0", new ReferenceInfoVO("18:0/18:0",28.9d,true,false));
    species.put("36:0", molSpecies);    
    molSpecies = new LinkedHashMap<String,ReferenceInfoVO>();
    molSpecies.put("18:0/18:1", new ReferenceInfoVO("18:0/18:1",27.2d,true,true));
    species.put("36:1", molSpecies);
    return species;
  }
  
  public static LinkedHashMap<String,LinkedHashMap<String,ReferenceInfoVO>> getDGSpeciesOrbitrap(){
    LinkedHashMap<String,LinkedHashMap<String,ReferenceInfoVO>> species = FoundBiologicalSpecies.getDGSpeciesOrbitrap();
    LinkedHashMap<String,ReferenceInfoVO> molSpecies = new LinkedHashMap<String,ReferenceInfoVO>();
    molSpecies.put("12:0/12:0/-", new ReferenceInfoVO("12:0/12:0/-",20.1d,true,false));
    species.put("24:0", molSpecies);
    molSpecies = new LinkedHashMap<String,ReferenceInfoVO>();
    molSpecies.put("20:0/18:0/-", new ReferenceInfoVO("20:0/18:0/-",36.1d,true,false));
    species.put("38:0", molSpecies);
    return species;
  }
  
  public static LinkedHashMap<String,LinkedHashMap<String,ReferenceInfoVO>> getTGSpeciesOrbitrap(){
    LinkedHashMap<String,LinkedHashMap<String,ReferenceInfoVO>> species = FoundBiologicalSpecies.getTGSpeciesOrbitrap();    
    return species;
  }
  
  public static LinkedHashMap<String,LinkedHashMap<String,ReferenceInfoVO>> getSMSpeciesOrbitrap(){
    LinkedHashMap<String,LinkedHashMap<String,ReferenceInfoVO>> species = new LinkedHashMap<String,LinkedHashMap<String,ReferenceInfoVO>>();
    LinkedHashMap<String,ReferenceInfoVO> molSpecies = new LinkedHashMap<String,ReferenceInfoVO>();   
    molSpecies = new LinkedHashMap<String,ReferenceInfoVO>();
    molSpecies.put("18:1;2/14:0", new ReferenceInfoVO("18:1;2/14:0",19.19d,true,false));
    species.put("32:1;2", molSpecies);
    molSpecies = new LinkedHashMap<String,ReferenceInfoVO>();
    molSpecies.put("18:1;2/15:0", new ReferenceInfoVO("18:1;2/15:0",20.67d,true,false));
    species.put("33:1;2", molSpecies);
    molSpecies = new LinkedHashMap<String,ReferenceInfoVO>();
    molSpecies.put("18:1;2/16:0", new ReferenceInfoVO("18:1;2/16:0",22.07d,true,false));
    species.put("34:1;2", molSpecies);
    molSpecies = new LinkedHashMap<String,ReferenceInfoVO>();
    molSpecies.put("18:1;2/16:1", new ReferenceInfoVO("18:1;2/16:1",19.98d,true,false));
    species.put("34:2;2", molSpecies);
    molSpecies = new LinkedHashMap<String,ReferenceInfoVO>();
    molSpecies.put("18:1;2/17:0", new ReferenceInfoVO("18:1;2/17:0",23.51d,true,false));
    species.put("35:1;2", molSpecies);
    molSpecies = new LinkedHashMap<String,ReferenceInfoVO>();
    molSpecies.put("18:1;2/18:0", new ReferenceInfoVO("18:1;2/18:0",24.82d,true,false));
    species.put("36:1;2", molSpecies);
    molSpecies = new LinkedHashMap<String,ReferenceInfoVO>();
    molSpecies.put("18:1;2/18:1", new ReferenceInfoVO("18:1;2/18:1",22.9d,true,false));
    species.put("36:2;2", molSpecies);
    molSpecies = new LinkedHashMap<String,ReferenceInfoVO>();
    molSpecies.put("18:1;2/20:0", new ReferenceInfoVO("18:1;2/20:0",27.46d,true,false));
    species.put("38:1;2", molSpecies);
    molSpecies = new LinkedHashMap<String,ReferenceInfoVO>();
    molSpecies.put("18:1;2/21:0", new ReferenceInfoVO("18:1;2/21:0",28.8d,true,false));
    species.put("39:1;2", molSpecies);
    molSpecies = new LinkedHashMap<String,ReferenceInfoVO>();
    molSpecies.put("18:1;2/22:0", new ReferenceInfoVO("18:1;2/22:0",29.87d,true,false));
    species.put("40:1;2", molSpecies);
    molSpecies = new LinkedHashMap<String,ReferenceInfoVO>();
    molSpecies.put("18:1;2/22:1", new ReferenceInfoVO("18:1;2/22:1",27.75d,true,false));
    species.put("40:1;2", molSpecies);
    molSpecies = new LinkedHashMap<String,ReferenceInfoVO>();
    molSpecies.put("18:1;2/23:0", new ReferenceInfoVO("18:1;2/23:0",31.01d,true,false));
    species.put("41:1;2", molSpecies);
    molSpecies = new LinkedHashMap<String,ReferenceInfoVO>();
    molSpecies.put("18:1;2/23:1", new ReferenceInfoVO("18:1;2/23:1",28.89d,true,false));
    species.put("41:2;2", molSpecies);
    molSpecies = new LinkedHashMap<String,ReferenceInfoVO>();
    molSpecies.put("18:1;2/24:0", new ReferenceInfoVO("18:1;2/24:0",32.06d,true,false));
    species.put("42:1;2", molSpecies);
    molSpecies = new LinkedHashMap<String,ReferenceInfoVO>();
    molSpecies.put("18:1;2/24:1", new ReferenceInfoVO("18:1;2/24:1",29.92d,true,false));
    species.put("42:2;2", molSpecies);
    molSpecies = new LinkedHashMap<String,ReferenceInfoVO>();
    molSpecies.put("18:1;2/24:2", new ReferenceInfoVO("18:1;2/24:2",28.19d,true,false));
    species.put("42:3;2", molSpecies);
    molSpecies = new LinkedHashMap<String,ReferenceInfoVO>();
    molSpecies.put("18:1;2/25:0", new ReferenceInfoVO("18:1;2/25:0",33.04d,true,false));
    species.put("43:1;2", molSpecies);
    molSpecies = new LinkedHashMap<String,ReferenceInfoVO>();
    molSpecies.put("18:1;2/25:1", new ReferenceInfoVO("18:1;2/25:1",30.98d,true,false));
    species.put("43:2;2", molSpecies);
    return species;
  }

  public static LinkedHashMap<String,LinkedHashMap<String,ReferenceInfoVO>> getCerSpeciesOrbitrap(){
    LinkedHashMap<String,LinkedHashMap<String,ReferenceInfoVO>> species = new LinkedHashMap<String,LinkedHashMap<String,ReferenceInfoVO>>();
    LinkedHashMap<String,ReferenceInfoVO> molSpecies = new LinkedHashMap<String,ReferenceInfoVO>();
    molSpecies.put("18:1;2/16:0", new ReferenceInfoVO("18:1;2/16:0",26.1d,true,false));
    species.put("34:1;2", molSpecies);
    molSpecies = new LinkedHashMap<String,ReferenceInfoVO>();
    molSpecies.put("18:1;2/17:0", new ReferenceInfoVO("18:1;2/17:0",27.3d,true,false));
    species.put("34:1;2", molSpecies);
    molSpecies = new LinkedHashMap<String,ReferenceInfoVO>();
    molSpecies.put("18:1;2/18:0", new ReferenceInfoVO("18:1;2/18:0",28.6d,true,false));
    species.put("36:1;2", molSpecies);
    molSpecies = new LinkedHashMap<String,ReferenceInfoVO>();
    molSpecies.put("18:1;2/18:1", new ReferenceInfoVO("18:1;2/18:1",26.7d,true,false));
    species.put("36:2;2", molSpecies);
    molSpecies = new LinkedHashMap<String,ReferenceInfoVO>();
    molSpecies.put("18:0;2/20:0", new ReferenceInfoVO("18:0;2/20:0",31.3d,true,false));
    species.put("38:0;2", molSpecies);
    molSpecies = new LinkedHashMap<String,ReferenceInfoVO>();
    molSpecies.put("18:1;2/20:0", new ReferenceInfoVO("18:1;2/20:0",30.8d,true,false));
    species.put("38:1;2", molSpecies);
    molSpecies = new LinkedHashMap<String,ReferenceInfoVO>();
    molSpecies.put("18:1;2/21:0", new ReferenceInfoVO("18:1;2/21:0",31.9d,true,false));
    species.put("39:1;2", molSpecies);
    molSpecies = new LinkedHashMap<String,ReferenceInfoVO>();
    molSpecies.put("18:1;2/22:0", new ReferenceInfoVO("18:1;2/22:0",32.9d,true,false));
    species.put("40:1;2", molSpecies);
    molSpecies = new LinkedHashMap<String,ReferenceInfoVO>();
    molSpecies.put("18:1;2/22:1", new ReferenceInfoVO("18:1;2/22:1",31.1d,true,false));
    species.put("40:2;2", molSpecies);
    molSpecies = new LinkedHashMap<String,ReferenceInfoVO>();
    molSpecies.put("18:1;2/23:0", new ReferenceInfoVO("18:1;2/23:0",33.8d,true,false));
    species.put("41:1;2", molSpecies);
    molSpecies = new LinkedHashMap<String,ReferenceInfoVO>();
    molSpecies.put("18:1;2/23:1", new ReferenceInfoVO("18:1;2/23:1",32.1d,true,false));
    species.put("41:2;2", molSpecies);
    molSpecies = new LinkedHashMap<String,ReferenceInfoVO>();
    molSpecies.put("18:1;2/24:0", new ReferenceInfoVO("18:1;2/24:0",34.4d,true,false));
    species.put("42:1;2", molSpecies);
    molSpecies = new LinkedHashMap<String,ReferenceInfoVO>();
    molSpecies.put("18:1;2/24:1", new ReferenceInfoVO("18:1;2/24:1",33.0d,true,false));
    species.put("42:2;2", molSpecies);
    molSpecies = new LinkedHashMap<String,ReferenceInfoVO>();
    molSpecies.put("18:1;2/24:2", new ReferenceInfoVO("18:1;2/24:2",31.5d,true,false));
    species.put("42:3;2", molSpecies);
    molSpecies = new LinkedHashMap<String,ReferenceInfoVO>();
    molSpecies.put("18:1;2/25:0", new ReferenceInfoVO("18:1;2/25:0",35.6d,true,false));
    species.put("43:1;2", molSpecies);
    molSpecies = new LinkedHashMap<String,ReferenceInfoVO>();
    molSpecies.put("18:1;2/25:1", new ReferenceInfoVO("18:1;2/25:1",33.9d,true,false));
    species.put("43:2;2", molSpecies);
    molSpecies = new LinkedHashMap<String,ReferenceInfoVO>();
    molSpecies.put("18:1;2/26:0", new ReferenceInfoVO("18:1;2/26:0",36.5d,true,false));
    species.put("44:1;2", molSpecies);
    return species;
  }
  
  
}
