package org.gwtvisualizationwrappers.client.biodalliance;

/*
 * #%L
 * GwtCytoscapeJs
 * %%
 * Copyright (C) 2015 GwtCytoscapeJs
 * %%
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *      http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 * #L%
 */
import com.google.gwt.json.client.JSONObject;
import com.google.gwt.json.client.JSONParser;
import com.google.gwt.json.client.JSONString;

/**
 * Source object represents a Biodalliance track source. 
 */
public class BiodallianceSource {
	String sourceName, sourceURI, entityId, indexEntityId, indexSourceURI;
	Long version, indexVersion;
	public static final String DEFAULT_STYLE_TYPE = "default";
	public static final String DEFAULT_STYLE_GLYPH_TYPE = "HISTOGRAM";
	public static final String DEFAULT_STYLE_COLOR = "#808080"; //grey
	public static final Integer DEFAULT_HEIGHT = 30;
	
	String styleType, styleGlyphType, styleColor;
	int trackHeightPx;

	public enum SourceType {
		BIGWIG, VCF, BED
	}
	SourceType sourceType;
	
	/**
	 * json keys
	 */
	public static final String SOURCE_NAME_KEY = "name";
	public static final String SOURCE_TYPE = "type";
	public static final String SOURCE_ENTITY_ID_KEY = "entityId";
	public static final String SOURCE_ENTITY_VERSION_KEY = "entityVersion";
	public static final String SOURCE_INDEX_ENTITY_ID_KEY = "indexEntityId";
	public static final String SOURCE_INDEX_ENTITY_VERSION_KEY = "indexEntityVersion";
	public static final String STYLE_TYPE_KEY = "styleType";
	public static final String STYLE_GLYPH_TYPE_KEY = "styleGlyphType";
	public static final String STYLE_COLOR_KEY = "color";
	public static final String STYLE_HEIGHT = "height";
	
	public BiodallianceSource() {
		initDefaults();
	}
	
	public BiodallianceSource(String json) {
		initializeFromJson(json);
	}
	
	
	public void initDefaults() {
		styleType = DEFAULT_STYLE_TYPE;
		styleGlyphType = DEFAULT_STYLE_GLYPH_TYPE;
		styleColor = DEFAULT_STYLE_COLOR;
		trackHeightPx = DEFAULT_HEIGHT;
	}
	
	public JSONObject toJsonObject() {
		JSONObject jsonObject = new JSONObject();
		jsonObject.put(SOURCE_NAME_KEY, new JSONString(sourceName));
		jsonObject.put(SOURCE_ENTITY_ID_KEY, new JSONString(entityId));
		if (version != null) {
			jsonObject.put(SOURCE_ENTITY_VERSION_KEY, new JSONString(version.toString()));	
		}
		if (indexEntityId != null) {
			jsonObject.put(SOURCE_INDEX_ENTITY_ID_KEY, new JSONString(indexEntityId));	
		}
		if (indexVersion != null) {
			jsonObject.put(SOURCE_INDEX_ENTITY_VERSION_KEY, new JSONString(indexVersion.toString()));	
		}
		jsonObject.put(STYLE_TYPE_KEY, new JSONString(styleType));
		jsonObject.put(STYLE_GLYPH_TYPE_KEY, new JSONString(styleGlyphType));
		jsonObject.put(STYLE_COLOR_KEY, new JSONString(styleColor));
		jsonObject.put(SOURCE_TYPE, new JSONString(sourceType.name()));
		jsonObject.put(STYLE_HEIGHT, new JSONString(Integer.toString(trackHeightPx)));
		return jsonObject;
	}
	
	public void initializeFromJson(String json) {
		JSONObject value = (JSONObject)JSONParser.parseStrict(json);
		setSourceName(value.get(SOURCE_NAME_KEY).isString().stringValue());
		String entityId = value.get(SOURCE_ENTITY_ID_KEY).isString().stringValue();
		Long version = null;
		if (value.containsKey(SOURCE_ENTITY_VERSION_KEY)) {
			String versionString = value.get(SOURCE_ENTITY_VERSION_KEY).isString().stringValue();
			version = Long.parseLong(versionString);
		}
		setEntity(entityId, version);
		
		String indexEntityId = null;
		if (value.containsKey(SOURCE_INDEX_ENTITY_ID_KEY)) {
			indexEntityId = value.get(SOURCE_INDEX_ENTITY_ID_KEY).isString().stringValue();
		}
		Long indexVersion = null;
		if (value.containsKey(SOURCE_INDEX_ENTITY_VERSION_KEY)) {
			String versionString = value.get(SOURCE_INDEX_ENTITY_VERSION_KEY).isString().stringValue();
			indexVersion = Long.parseLong(versionString);
		}
		
		setIndexEntity(indexEntityId, indexVersion);
		
		setStyleType(value.get(STYLE_TYPE_KEY).isString().stringValue());
		setStyleGlyphType(value.get(STYLE_GLYPH_TYPE_KEY).isString().stringValue());
		setStyleColor(value.get(STYLE_COLOR_KEY).isString().stringValue());
		setTrackHeightPx(Integer.parseInt(value.get(STYLE_HEIGHT).isString().stringValue()));
		String sourceTypeString = value.get(SOURCE_TYPE).isString().stringValue();
		setSourceType(SourceType.valueOf(sourceTypeString));
	}
	
	public void setStyleType(String styleType) {
		this.styleType = styleType;
	}

	public void setStyleGlyphType(String styleGlyphType) {
		this.styleGlyphType = styleGlyphType;
	}

	public void setStyleColor(String styleColor) {
		this.styleColor = styleColor;
	}
	
	public String getSourceName() {
		return sourceName;
	}
	public String getSourceURI() {
		return sourceURI;
	}
	public String getSourceIndexURI() {
		return indexSourceURI;
	}
	public String getStyleType() {
		return styleType;
	}
	public String getStyleGlyphType() {
		return styleGlyphType;
	}
	public String getStyleColor() {
		return styleColor;
	}
	public int getTrackHeightPx() {
		return trackHeightPx;
	}
	
	public SourceType getSourceType() {
		return sourceType;
	}
	
	public String getEntityId() {
		return entityId;
	}

	public void setEntity(String entityId, Long version) {
		this.entityId = entityId;
		this.version = version;
	}

	public String getIndexEntityId() {
		return indexEntityId;
	}

	public void setIndexEntity(String indexEntityId, Long indexVersion) {
		this.indexEntityId = indexEntityId;
		this.indexVersion = indexVersion;
	}
	
	public Long getVersion() {
		return version;
	}
	
	public Long getIndexVersion() {
		return indexVersion;
	}

	public void setSourceName(String sourceName) {
		this.sourceName = sourceName;
	}

	public void setTrackHeightPx(int trackHeightPx) {
		this.trackHeightPx = trackHeightPx;
	}
	
	public void setSourceType(SourceType sourceType) {
		this.sourceType = sourceType;
	}
	
	public void setSourceURI(String sourceURI) {
		this.sourceURI = sourceURI;
	}
	
	public void setIndexSourceURI(String indexSourceURI) {
		this.indexSourceURI = indexSourceURI;
	}
}
