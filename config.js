// RethinkDB settings
exports.host = process.env.RETHINK_HOST;	   					// RethinkDB host
exports.authKey = process.env.RETHINK_AUTH_KEY;   				// RethinkDB authentification key
exports.port = 28015;          									// RethinkDB driver port

// Express settings
exports.expressPort = 3000;    // Port used by express
exports.debug = true;          // Debug mode
exports.network = '127.0.0.1'  // Network the node app will run on
