// ***************************//
// Menu manipulation section
// ***************************//

/**
 * Function to set the iframe with the network HTML
 * 
 * @param {*} networkCode_ 
 * @param {*} networkIndex_ 
 */
function setNetwork(networkCode_, networkIndex_) {
    // Set the HTML to the iframe
    document.getElementById('networkIframe').src = 'https://igorabrandao.com.br/kegg-pathway-bottleneck/networks/' +
        networkIndex_ + '_' + networkCode_ + '.html';
}

function setNetworkAttributes(menuElement_, title_) {
    // Set the network title in HTML
    document.getElementById('graphTitle').innerHTML = 'Network dynamic visualization';
    //document.getElementById('graphTitle').innerHTML = title_;

    // Remove the highlight from all menus
    var elems = document.querySelectorAll(".highlight");

    [].forEach.call(elems, function(el) {
        el.classList.remove("highlight");
    });

    // Add the highlight class to the selected menu
    document.getElementById(menuElement_).className += " highlight";
}

// *************************//
// Loader handler section
// *************************//

/**
 * Check if the page is ready
 * 
 * @param {*} callback 
 */
function onReady(callback) {
    var intervalID = window.setInterval(checkReady, 1000);

    function checkReady() {
        if (document.getElementsByTagName('body')[0] !== undefined) {
            window.clearInterval(intervalID);
            callback.call(this);
        }
    }
}

/**
 * Show the loader
 * 
 * @param {*} id 
 * @param {*} value 
 */
function show(id, value) {
    document.getElementById(id).style.display = value ? 'block' : 'none';
}

/**
 * Function to hide the loader
 */
onReady(function () {
    show('page', true);
    show('loading', false);
});