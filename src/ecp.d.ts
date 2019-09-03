import {BIG} from "./big";
import {ECP2} from "./ecp2";
import {EcpCommon} from "./ecp-common";
import {FP} from "./fp";

export class ECP extends EcpCommon<ECP> {

    x: FP;
    y: FP;
    z: FP;

    public constructor(input?: ECP);

    public setxy(x: BIG, y: BIG): void;
    public setx(x: BIG): void;
    public copy(p: ECP): void;
    public neg(): void;
    public inf(): void;
    public affine(): void;
    public getX(): BIG;
    public getY(): BIG;
    public getS(): boolean;
    public sub(p: ECP|ECP2): void;
    public getX(): ECP;
    public toBytes(destination: Uint8Array, compress: boolean): void;

    static generator(): ECP;
    static fromBytes(array: Uint8Array): ECP;
}
