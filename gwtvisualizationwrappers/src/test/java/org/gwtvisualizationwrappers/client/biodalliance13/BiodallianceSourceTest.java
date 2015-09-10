package org.gwtvisualizationwrappers.client.biodalliance13;

import static org.junit.Assert.*;

import org.junit.Ignore;
import org.junit.Test;

import com.google.gwt.json.client.JSONObject;

@Ignore
public class BiodallianceSourceTest  {

	@Test
	public void testInitDefaults() {
		BiodallianceSource source = new BiodallianceSource();
		assertEquals(BiodallianceSource.DEFAULT_STYLE_TYPE, source.getStyleType());
	}

	@Test
	public void testJsonRoundTrip() {
		//create a source, write it out as json, initialize another source with that json,
		// and then verify that the two objects are equal.
		BiodallianceSource s1 = new BiodallianceSource();
		s1.setEntity("syn22", 2L);
		s1.setHeightPx(20);
		s1.setIndexEntity("syn44", 3L);
		s1.setSourceName("a round trip json test");
		s1.setStyleColor("red");
		s1.setStyleGlyphType("PLIMSOLL");
		s1.setStyleType("default");
		
		JSONObject jsonObject = s1.toJsonObject();
		
		BiodallianceSource s2 = new BiodallianceSource();
		s2.initializeFromJson(jsonObject.toString());
		assertEquals(s1, s2);
	}

}
