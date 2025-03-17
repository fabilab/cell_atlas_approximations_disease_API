const atlasapprox_disease = require("../index");

const delay = (ms) => new Promise(resolve => setTimeout(resolve, ms));

describe('Differential tests', () => {
    beforeEach(async () => {
        // Wait for 1s before each test
        await delay(1000);
    });

    test("differential gene expression returns expected structure", async () => {
        const result = await atlasapprox_disease.differential_gene_expression({
            tissue: "eye",
        });
        expect(result).toBeDefined();
        expect(Array.isArray(result)).toBe(true);
        result.forEach(entry => {
            expect(entry).toHaveProperty("gene");
            expect(entry).toHaveProperty("state_expr");
            expect(entry).toHaveProperty("baseline_expr");
            expect(entry).toHaveProperty("state_fraction");
            expect(entry).toHaveProperty("baseline_fraction");
            expect(typeof entry.gene).toBe("string");
            expect(typeof entry.state_expr).toBe("number");
            expect(typeof entry.baseline_expr).toBe("number");
            expect(typeof entry.state_fraction).toBe("number");
            expect(typeof entry.baseline_fraction).toBe("number");
            expect(entry.baseline.toLowerCase()).toBe("normal");
        });
    }, 50000);

    test("differential gene expression with multiple filters", async () => {
        const result = await atlasapprox_disease.differential_gene_expression({
            tissue: "pancreas",
            differential_axis: "sex",
            cell_type: "endothelial"
        });
        expect(result).toBeDefined();
        expect(Array.isArray(result)).toBe(true);
        result.forEach(entry => {
            expect(entry.differential_axis.toLowerCase().includes("sex")).toBeTruthy();
            expect(entry.cell_type.toLowerCase().includes("endothelial")).toBeTruthy();
        });
    }, 50000);

    test("differential cell type abundance returns expected structure", async () => {
        const result = await atlasapprox_disease.differential_cell_type_abundance({
            tissue: "eye",
            disease: "cataract"
        });
        expect(result).toBeDefined();
        expect(Array.isArray(result)).toBe(true);
        result.forEach(entry => {
            expect(entry).toHaveProperty("cell_type");
            expect(entry).toHaveProperty("frac_baseline");
            expect(entry).toHaveProperty("frac_disease");
            expect(entry).toHaveProperty("delta_frac");
            expect(typeof entry.cell_type).toBe("string");
            expect(typeof entry.frac_baseline).toBe("number");
            expect(typeof entry.frac_disease).toBe("number");
            expect(typeof entry.delta_frac).toBe("number");
            expect(entry.baseline.toLowerCase()).toBe("normal");
        });
    }, 50000);

    test("differential cell type with filters", async () => {
        const result = await atlasapprox_disease.differential_cell_type_abundance({
            tissue: "lung",
            cell_type: "macrophage",
            sex: "female",
        });

        expect(result).toBeDefined();
        expect(Array.isArray(result)).toBe(true);
        result.forEach(entry => {
            expect(entry.tissue_general.toLowerCase().includes("lung")).toBeTruthy();
            expect(entry.cell_type.toLowerCase().includes("macrophage")).toBeTruthy();
            expect(entry.sex.toLowerCase().includes("female")).toBeTruthy();
        });
    }, 50000);
});
