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
function createPlainObject(cytoscapeGraphJson, containerElement, readyFunction){
	var obj = JSON.parse(cytoscapeGraphJson);
	obj.container=containerElement;
	obj.ready=readyFunction;
	
	return obj;
};

