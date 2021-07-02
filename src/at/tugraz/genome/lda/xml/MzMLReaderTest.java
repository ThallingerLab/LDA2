package at.tugraz.genome.lda.xml;

import static org.junit.jupiter.api.Assertions.*;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import java.lang.reflect.Field;

class MzMLReaderTest
{

  @Test
  @DisplayName("Sets a new AddScan interface correctly.")
  void setAddersTest() throws NoSuchFieldException, IllegalAccessException
  {
    //given
    AddScan[] adder = new AddScan[1];
    final XmlSpectraReader reader = new MzMLReader(adder, false);
    
    //when
    AddScan[] newAdder = new AddScan[1];
    reader.setAdders(newAdder);
    
    //then
    final Field field = reader.getClass().getSuperclass().getDeclaredField("adders_");
    field.setAccessible(true);
    assertEquals(newAdder, field.get(reader));
  }
  
  @Test
  @DisplayName("Gets parseMsMs_ correctly.")
  void getParseMsMsTest() throws NoSuchFieldException, IllegalAccessException
  {
    //given
    AddScan[] adder = new AddScan[1];
    final MzMLReader reader = new MzMLReader(adder, false);
    
    final Field field = reader.getClass().getSuperclass().getDeclaredField("parseMsMs_");
    field.setAccessible(true);
    field.set(reader, true);

    //when
    final boolean result = reader.getParseMsMs();

    //then
    assertEquals(true, result);
  }

  @Test
  @DisplayName("Gets multiplicationFactorForInt_ correctly.")
  void getMultiplicationFactorForIntTest() throws NoSuchFieldException, IllegalAccessException
  {
    //given
    AddScan[] adder = new AddScan[1];
    final MzMLReader reader = new MzMLReader(adder, false);
    
    final Field field = reader.getClass().getSuperclass().getDeclaredField("multiplicationFactorForInt_");
    field.setAccessible(true);
    field.set(reader, 8);

    //when
    final int result = reader.getMultiplicationFactorForInt();

    //then
    assertEquals(8, result);
  }
  
  @Test
  @DisplayName("Sets highestMz_ correctly.")
  void setHighestMzTest() throws NoSuchFieldException, IllegalAccessException
  {
    //given
    AddScan[] adder = new AddScan[1];
    final MzMLReader reader = new MzMLReader(adder, false);
    
    //when
    reader.setHighestMz(3);
    
    //then
    final Field field = reader.getClass().getSuperclass().getDeclaredField("highestMz_");
    field.setAccessible(true);
    assertEquals(3, field.get(reader));
  }
  
  @Test
  @DisplayName("Gets highestMz_ correctly.")
  void getHighestMzTest() throws NoSuchFieldException, IllegalAccessException
  {
    //given
    AddScan[] adder = new AddScan[1];
    final MzMLReader reader = new MzMLReader(adder, false);
    
    final Field field = reader.getClass().getSuperclass().getDeclaredField("highestMz_");
    field.setAccessible(true);
    field.set(reader, 3);

    //when
    final int result = reader.getHighestMz();

    //then
    assertEquals(3, result);
  }
  
  @Test
  @DisplayName("Sets lowestMz_ correctly.")
  void setLowestMzTest() throws NoSuchFieldException, IllegalAccessException
  {
    //given
    AddScan[] adder = new AddScan[1];
    final MzMLReader reader = new MzMLReader(adder, false);
    
    //when
    reader.setLowestMz(5);
    
    //then
    final Field field = reader.getClass().getSuperclass().getDeclaredField("lowestMz_");
    field.setAccessible(true);
    assertEquals(5, field.get(reader));
  }
  
  @Test
  @DisplayName("Gets lowestMz_ correctly.")
  void getLowestMzTest() throws NoSuchFieldException, IllegalAccessException
  {
    //given
    AddScan[] adder = new AddScan[1];
    final MzMLReader reader = new MzMLReader(adder, false);
    
    final Field field = reader.getClass().getSuperclass().getDeclaredField("lowestMz_");
    field.setAccessible(true);
    field.set(reader, 1);

    //when
    final int result = reader.getLowestMz();

    //then
    assertEquals(1, result);
  }
  
  @Test
  @DisplayName("Sets polaritySwitching_ correctly.")
  void setPolaritySwitchingTest() throws NoSuchFieldException, IllegalAccessException
  {
    //given
    AddScan[] adder = new AddScan[1];
    final MzMLReader reader = new MzMLReader(adder, false);
    
    //when
    reader.setPolaritySwitching(true);
    
    //then
    final Field field = reader.getClass().getSuperclass().getDeclaredField("polaritySwitching_");
    field.setAccessible(true);
    assertEquals(true, field.get(reader));
  }
  
  @Test
  @DisplayName("Gets polaritySwitching_ correctly.")
  void usesPolaritySwitchingTest() throws NoSuchFieldException, IllegalAccessException
  {
    //given
    AddScan[] adder = new AddScan[1];
    final XmlSpectraReader reader = new MzMLReader(adder, false);
    
    final Field field = reader.getClass().getSuperclass().getDeclaredField("polaritySwitching_");
    field.setAccessible(true);
    field.set(reader, true);

    //when
    final boolean result = reader.usesPolaritySwitching();

    //then
    assertEquals(true, result);
  }
  

  
  
}
