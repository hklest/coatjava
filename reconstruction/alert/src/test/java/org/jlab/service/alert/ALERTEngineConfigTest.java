package org.jlab.service.alert;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.HashMap;
import java.util.Map;
import org.jlab.rec.ahdc.GlobalFit.GFitConfig;
import org.junit.Test;

/** YAML keys of the ALERT global track fit. */
public class ALERTEngineConfigTest {

  private static ALERTEngine configured(Map<String, String> config, boolean[] ok) {
    ALERTEngine engine = new ALERTEngine();
    ok[0] = engine.readGfitConfig(config::get);
    return engine;
  }

  @Test
  public void gfitOffByDefault() {
    boolean[] ok = new boolean[1];
    ALERTEngine engine = configured(new HashMap<>(), ok);
    assertTrue(ok[0]);
    assertFalse(engine.isGfitOn());
    assertTrue(engine.isGfitRecovery());
    assertFalse(engine.isGfitTime());
    assertFalse(engine.isGfitDedx());
  }

  @Test
  public void gfitSwitchesParsed() {
    Map<String, String> config = new HashMap<>();
    config.put("gfit", "true");
    config.put("gfitRecovery", "false");
    config.put("gfitTime", "true");
    config.put("gfitDedx", "true");
    boolean[] ok = new boolean[1];
    ALERTEngine engine = configured(config, ok);
    assertTrue(ok[0]);
    assertTrue(engine.isGfitOn());
    assertFalse(engine.isGfitRecovery());
    assertTrue(engine.isGfitTime());
    assertTrue(engine.isGfitDedx());
  }

  @Test
  public void gfitBadValue() {
    for (String key : new String[] {"gfit", "gfitRecovery", "gfitTime", "gfitDedx"}) {
      Map<String, String> config = new HashMap<>();
      config.put(key, "maybe");
      boolean[] ok = new boolean[1];
      configured(config, ok);
      assertFalse(key, ok[0]);
    }
  }

  @Test
  public void gfitSearchValues() {
    boolean[] ok = new boolean[1];
    assertEquals(GFitConfig.Search.FAST, configured(new HashMap<>(), ok).getGfitSearch());
    assertTrue(ok[0]);
    String[][] cases = {{"fast", "FAST"}, {"full", "FULL"}};
    for (String[] c : cases) {
      Map<String, String> config = new HashMap<>();
      config.put("gfitSearch", c[0]);
      ALERTEngine engine = configured(config, ok);
      assertTrue(c[0], ok[0]);
      assertEquals(GFitConfig.Search.valueOf(c[1]), engine.getGfitSearch());
    }
    Map<String, String> config = new HashMap<>();
    config.put("gfitSearch", "quick");
    configured(config, ok);
    assertFalse(ok[0]);
  }
}
