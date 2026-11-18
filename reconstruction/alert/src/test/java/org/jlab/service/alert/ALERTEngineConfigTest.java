package org.jlab.service.alert;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.HashMap;
import java.util.Map;
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
}
