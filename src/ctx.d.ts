import {BIG, DBIG} from "./big";
import {FP} from "./fp";
import {FP2} from "./fp2";
import {FP4} from "./fp4";
import {FP12} from "./fp12";
import {ECP} from "./ecp";
import {ECP2} from "./ecp2";
import {PAIR} from "./pair";
import {RomField} from "./rom_field";

export class CTX {
    public BIG: typeof BIG;
    public DBIG: typeof DBIG;
    public FP: typeof FP;
    public FP2: typeof FP2;
    public FP4: typeof FP4;
    public FP12: typeof FP12;
    public ECP: typeof ECP;
    public ECP2: typeof ECP2;

    public PAIR: PAIR;
    public ROM_FIELD: RomField;
}
