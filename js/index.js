const api_version = "v1";
const api_uri_default = "https://api-disease.atlasapprox.org/" + api_version + "/";
// const api_uri_default = "http://127.0.0.1:5000/" + api_version + "/";  // Local API
let api_uri = api_uri_default;

const setAPIURI = (uri) => {
  api_uri = uri;
};

const resetAPIURI = () => {
  api_uri = api_uri_default;
};

async function _callEndpoint(endpoint, params = {}, method = "GET") {
  Object.keys(params).forEach(key => !params[key] && delete params[key]);

  let uri = api_uri + endpoint;
  let options = {
      method: method,
      headers: {
          "Content-Type": "application/json"
      }
  };

  const uriSuffix = new URLSearchParams(params).toString();
  if (uriSuffix !== "") uri += "?" + uriSuffix;

  let response = await fetch(uri, options);

  const data = await response.json();

  if (!response.ok) {
      throw {
          status: response.status,
          message: data.message,
          error: data.error,
      };
  }
  return data;
}

let isString = value => typeof value === 'string' || value instanceof String;

const formatFeatures = (features) => {
  if ((features === undefined) || (features === null))
    return features

  if (!isString(features))
    return features.join(",");

  return features;
};

///////////////////////////////////////////////////////////////////////////////
// Public API Functions
///////////////////////////////////////////////////////////////////////////////
async function metadata({ disease = null, cell_type = null, tissue = null, sex = null, development_stage = null }) {
  let params = { disease, cell_type, tissue, sex, development_stage };
  return await _callEndpoint("metadata", params=params);
}

async function differential_cell_type_abundance({ differential_axis = "disease", disease = null, cell_type = null, tissue = null, sex = null, development_stage = null }) {
  let params = { differential_axis, disease, cell_type, tissue, sex, development_stage }
  return await _callEndpoint("differential_cell_type_abundance", params=params, "POST");
}

async function differential_gene_expression({ differential_axis = "disease", disease = null, cell_type = null, tissue = null, sex = null, development_stage = null, top_n = null, feature = null, method = "delta_fraction" }) {
  let params = { differential_axis, disease, cell_type, tissue, sex, development_stage, top_n, feature, method }
  return await _callEndpoint("differential_gene_expression", params=params, "POST");
}

async function highest_measurement({ feature, number = null }) {
  let params = { feature, number }
  return await _callEndpoint("highest_measurement", params=params, "POST");
}

async function average({ features, disease = null, cell_type = null, tissue = null, sex = null, development_stage = null, unique_ids = null, include_normal = null }) {
  features = formatFeatures(features);
  let params = { features, disease, cell_type, tissue, sex, development_stage, unique_ids, include_normal };
  return await _callEndpoint("average", params=params, "POST");
}

async function fraction_detected({ features, disease = null, cell_type = null, tissue = null, sex = null, development_stage = null, unique_ids = null, include_normal = null }) {
  features = formatFeatures(features);
  let params = { features, disease, cell_type, tissue, sex, development_stage, unique_ids, include_normal };
  return await _callEndpoint("fraction_detected", params=params, "POST");
}

async function dotplot({ features, disease = null, cell_type = null, tissue = null, sex = null, development_stage = null, unique_ids = null, include_normal = null }) {
  features = formatFeatures(features);
  let params = { features, disease, cell_type, tissue, sex, development_stage, unique_ids, include_normal }
  return await _callEndpoint("dotplot", params, "POST");
}

///////////////////////////////////////////////////////////////////////////////
// Exported API Object
///////////////////////////////////////////////////////////////////////////////
const atlasapprox_disease = {
  metadata,
  differential_cell_type_abundance,
  differential_gene_expression,
  highest_measurement,
  average,
  fraction_detected,
  dotplot,
  setAPIURI,
  resetAPIURI,
  api_version,
};

module.exports = atlasapprox_disease;
