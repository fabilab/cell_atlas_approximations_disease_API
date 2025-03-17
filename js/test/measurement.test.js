const atlasapprox_disease = require("../index");

const delay = (ms) => new Promise(resolve => setTimeout(resolve, ms));

describe('Differential tests', () => {
    beforeEach(async () => {
        // Wait for 1s before each test
        await delay(1000);
    });
    
    test("average with multiple filters", async () => {
        const result = await atlasapprox_disease.average({
            features: "IL6, AGT",
            disease: "covid",
            cell_type: "fibroblast",
        });
        expect(result).toBeDefined();
        expect(Array.isArray(result)).toBe(true);
        result.forEach(entry => {
            // Ensure the disease field contains "covid"
            expect(entry.disease.toLowerCase().includes("covid")).toBeTruthy();

            expect(entry.cell_type.toLowerCase().includes("fibroblast")).toBeTruthy();

            // Ensure gene names exist as properties inside each entry
            ["IL6", "AGT"].forEach(gene => {
                expect(entry).toHaveProperty(gene);
                expect(typeof entry[gene]).toBe("number");
            });
        });
    }, 10000);

    test("fraction with multiple filters", async () => {
        const result = await atlasapprox_disease.average({
            features: "ANKRD12, SPP1",
            disease: "acute kidney failure",
        });
        expect(result).toBeDefined();
        expect(Array.isArray(result)).toBe(true);
        result.forEach(entry => {
            expect(entry.disease.toLowerCase().includes("acute kidney failure")).toBeTruthy();
            ["ANKRD12", "SPP1"].forEach(gene => {
                expect(entry).toHaveProperty(gene);
                expect(typeof entry[gene]).toBe("number");
            });
        });
    }, 10000);

    test("dotplot with multiple filters", async () => {
        const result = await atlasapprox_disease.dotplot({
            features: "ANKRD12, SPP1",
            disease: "acute kidney failure",
        });
        expect(result).toBeDefined();
        expect(Array.isArray(result)).toBe(true);
        result.forEach(entry => {
            expect(entry.disease.toLowerCase().includes("acute kidney failure")).toBeTruthy();
            ["ANKRD12", "SPP1"].forEach(gene => {
                expect(entry).toHaveProperty(gene);
                expect(entry).toHaveProperty(`fraction_${gene}`);
                expect(typeof entry[gene]).toBe("number");
            });
        });
    }, 10000);

    test("dotplot include normal", async () => {
        const result = await atlasapprox_disease.dotplot({
            features: "ANKRD12, SPP1",
            disease: "acute kidney failure",
            include_normal: true,
        });
        expect(result).toBeDefined();
        expect(Array.isArray(result)).toBe(true);
        
        result.forEach(entry => {
            expect(["acute kidney failure", "normal"].includes(entry.disease)).toBeTruthy();
            ["ANKRD12", "SPP1"].forEach(gene => {
                expect(entry).toHaveProperty(gene);
                expect(entry).toHaveProperty(`fraction_${gene}`);
                expect(typeof entry[gene]).toBe("number");
            });
        });
    }, 10000);

});