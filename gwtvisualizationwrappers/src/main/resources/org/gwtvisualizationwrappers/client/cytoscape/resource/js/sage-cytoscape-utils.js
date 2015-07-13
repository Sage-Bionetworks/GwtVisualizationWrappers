/**
 * Creating a js Object obj inside of a GWT JSNI method will cause obj.constructor === Object to evaluation to false.
 * Creating the js Object outside of that environment forces it to create a vanilla js Object.
 * This is needed for Cytoscape.js because it checks the input options Object in this way. 
 * @param containerElement
 * @param elementsArray
 * @param styleDefinition
 * @param readyFunction
 * @returns Object
 */
function createPlainObject(cyjs, styleJson, containerElement, readyFunction){
	var obj = JSON.parse(cyjs);
	obj.container=containerElement;
	obj.ready=readyFunction;
	
	try {
		var styles = JSON.parse(styleJson);
		//try to use the last style defined
		obj.style=styles[styles.length - 1].style;
	}
	catch(err) {
		console.log('error loading Cytoscape style property.', err);
	}
	return obj;
};

