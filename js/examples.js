const atlasapprox_disease = require("./index");

// Test cases for each function

async function testMetadata() {
    try {
        console.log("Testing /metadata...");
        const response = await atlasapprox_disease.metadata({ disease: "diabetes" });
        console.log("Metadata Response:", response.slice(0, 10));
    } catch (error) {
        console.error("Error:", error);
    }
}

async function testDifferentialCellTypeAbundance() {
    try {
        console.log("Testing /differential_cell_type_abundance...");
        const response = await atlasapprox_disease.differential_cell_type_abundance({
            disease: "can",
            // cell_type: "T cell",
        });
        console.log("Differential Cell Type Abundance Response:", response);
    } catch (error) {
        console.error("Error:", error);
    }
}

async function testDifferentialGeneExpression() {
    try {
        console.log("Testing /differential_gene_expression...");
        const response = await atlasapprox_disease.differential_gene_expression({
            disease: "Covid",
            cell_type: "T cell",
            top_n: 10
        });
        console.log("Differential Gene Expression Response:", response);
    } catch (error) {
        console.error("Error:", error);
    }
}

async function testHighestMeasurement() {
    try {
        console.log("Testing /highest_measurement...");
        const response = await atlasapprox_disease.highest_measurement({ feature: "IL6" });
        console.log("Highest Measurement Response:", response);
    } catch (error) {
        console.error("Error:", error);
    }
}

async function testAverage() {
    try {
        console.log("Testing /average...");
        const response = await atlasapprox_disease.average({
            features: "IL6, AGT",
            disease: "Covid",
            cell_type: "fibroblast",
            include_normal: true,
        });
        console.log("Average Expression Response:", response);
    } catch (error) {
        console.error("Error:", error);
    }
}

async function testAverageUniqueId() {
    try {
        console.log("Testing /average with unique ids...");
        const response = await atlasapprox_disease.average({
            features: "CD19, CD68",
            unique_ids: "f9ac955738b4d19f4565c3eebd38a1b6,"
        });
        console.log("Average Expression Response:", response);
    } catch (error) {
        console.error("Error:", error);
    }
}

async function testFractionDetected() {
    try {
        console.log("Testing /fraction_detected...");
        const response = await atlasapprox_disease.fraction_detected({
            features: "IL6, AGT",
            disease: "Covid",
            cell_type: "fibroblast",
            include_normal: true,
        });
        console.log("Fraction Detected Response:", response);
    } catch (error) {
        console.error("Error:", error);
    }
}

async function testDotplot() {
    try {
        console.log("Testing /dotplot...");
        const response = await atlasapprox_disease.dotplot({
            features: "IL6, AGT",
            disease: "Covid",
            cell_type: "fibroblast",
            include_normal: true,
        });
        console.log("Dotplot Response:", response);
    } catch (error) {
        console.error("Error:", error);
    }
}

// Run all test cases
(async () => {
    await testMetadata();
    // await testDifferentialCellTypeAbundance();
    // await testDifferentialGeneExpression();
    // await testHighestMeasurement();
    // await testAverage();
    await testAverageUniqueId();
    // await testFractionDetected();
    // await testDotplot();
})();
