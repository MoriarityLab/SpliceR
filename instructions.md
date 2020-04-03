## What this is

The purpose of this software is to predict base editing sgRNAs to use in the disruption of splice sites for gene modulation, or knockout. This approach is especially amenable to multiplex gene editing relative to Cas9 nuclease based approaches.

### Citation

If our software helps you -- please cite us!

Mitchell G. Kluesner†, Walker S. Lahr†, Cara-Lin Lonetree, Branden A. Smeester, Patricia N. Claudio-Vázquez, Samuel Pitzen, Madison J. Vignes, Samantha C. Lee, Samuel P. Bingea, Aneesha A. Andrews, Beau R. Webber†, Branden S. Moriarity†
**CRISPR-Cas9 cytidine and adenine base editing of splice-sites mediates highly-efficient disruption of proteins in primary cells**. *In review*.

† These authors contributed equally.

## Running the prediction software

To run this software, please provide:

*   An Ensembl transcript ID for a desired transcript to target
*   The PAM specificity of your desired base editor
*	The species being targeted
*	The desired type base editor
*	The desired splice-site to target

### Outputs

**Predicted guides**: This is a table of predicted sgRNAs based on the input parameters. Base editing efficiencies are predictions.

**Ensembl gene map**: This is a rendering of the Ensembl gene page for the input transcript ID. This table can be used to compare different isoforms of the target gene. 

**Download button**: This button located on the left sidebar panel can be used to download the table provided in the "Predicted guides" tab.
